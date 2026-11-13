#####################################################
#   Code to extract observation level data from LOTUS
#   Requires approved NDA access 
#   This version uses Free TDS as the driver
#
######################################################
library(DBI)
library(odbc)
library(glue)

driver.name <- "/opt/homebrew/opt/freetds/lib/libtdsodbc.so"
db.name<-"LOTUS"
port<-"11433"
user.name<-Sys.getenv("LODS_USER")#Notethedoublebackslashformac
pwd<-Sys.getenv("LODS_PASS") 
host.name<-Sys.getenv("LODS_LOC")

con <- dbConnect(
  odbc::odbc(),
  Driver = driver.name,
  Server = host.name,
  Port = port,
  Database = db.name,
  UID = user.name,
  PWD = pwd,
  TDS_Version = "auto",
  Protocol = "TCPIP",
  Encrypt = "yes",
  TrustServerCertificate = "no",
  ConnectionTimeout = 30
)

setwd("set your working directory here")

syr <- 2005
eyr <- 2025

start_date <- paste0(syr, "-01-01")
end_date   <- paste0(eyr + 1, "-01-01")

qtext <- glue_sql("
    SELECT
        YEAR(se.HAUL_BEGIN_DATETIME) AS year,
        tr.TRIP_NUM AS trip,

        COUNT(
            CASE
                WHEN ca.species_code = '506'
                THEN ca.ln_num
            END
        ) AS leatherback,

        COUNT(
            CASE
                WHEN ca.species_code IN ('504', '500')
                THEN ca.ln_num
            END
        ) AS loggerhead

    FROM LOTUS.newobs.trips tr

    INNER JOIN LOTUS.newobs.lop_declared_trip_type dtt
        ON tr.declared_trip_type_id = dtt.ID

    INNER JOIN LOTUS.newobs.lop_set_environ_v se
        ON tr.trip_num = se.t_trip_num

    LEFT JOIN LOTUS.newobs.catch ca
        ON se.t_trip_num = ca.t_trip_num
        AND se.set_num = ca.s_set_num

    WHERE
        tr.trip_num LIKE 'LL%'
        AND dtt.code = 'S'
        AND se.HAUL_BEGIN_DATETIME >= {start_date}
        AND se.HAUL_BEGIN_DATETIME < {end_date}

    GROUP BY
        YEAR(se.HAUL_BEGIN_DATETIME),
        tr.TRIP_NUM

    ORDER BY
        YEAR(se.HAUL_BEGIN_DATETIME),
        tr.TRIP_NUM
", .con = con)

dfile <- dbGetQuery(con, qtext)

colnames(dfile) <- c("year", "trip", "leatherback", "loggerhead")

# ---- 1. Trip-level aggregation ----
trip_turtles <- dfile %>%
  group_by(year, trip) %>%
  summarise(
    n_sets      = n(),
    leatherback = sum(leatherback, na.rm = TRUE),
    loggerhead  = sum(loggerhead, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(leatherback, loggerhead), names_to = "species", values_to = "count") %>%
  mutate(block = case_when(
    year >= 2005 & year <= 2020 ~ "Pre",
    year >= 2021 & year <= 2025 ~ "Post",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(block))

# ---- 2. Trip limits per species, and censoring flag ----
trip_limits <- tibble(
  species = c("leatherback", "loggerhead"),
  limit   = c(2, 5)
)

trip_turtles_cens <- trip_turtles %>%
  left_join(trip_limits, by = "species") %>%
  mutate(censored = if_else(block == "Post" & count >= limit,
                             "right", "none"))

trip_turtles_cens <- trip_turtles_cens %>%
  mutate(count_fit = if_else(censored == "right", count - 1, count)) # shift given the brms implentation of the NB> or >=

trip_turtles_cens <- trip_turtles_cens %>%
  mutate(block = factor(block, levels = c("Pre", "Post")))


save(trip_turtles_cens,file=f"trip_turtles_cens.Rdata")
