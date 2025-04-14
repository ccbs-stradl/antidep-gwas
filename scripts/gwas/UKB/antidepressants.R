library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(xml2)

# connect to database
con <- DBI::dbConnect(RSQLite::SQLite(), here::here('data/phenotypes/health/health_outcomes/healthoutcomes.db'))

# get registrations
gp_registrations <- tbl(con, 'gp_registrations')

# get prescribing records
gp_scripts <- tbl(con, 'gp_scripts')

# participants with records
eids_with_records <- gp_registrations %>% distinct(eid) %>% collect()

# antidepressants (BNF codes 04.04.01-04
gp_scripts_ad <- gp_scripts %>%
    filter(bnf_code %LIKE% '04.03.%') %>%
    select(eid, bnf_code, drug_name) %>%
    distinct() %>%
    collect()

# binary prescribed status
gp_prescribed_ad <-
gp_scripts_ad %>%
mutate(bnf=str_replace_all(str_sub(bnf_code, 1, 8), pattern="\\.", replacement="")) %>%
mutate(bnf=recode(bnf, "040301"="bnf040301_tca",
                        "040302"="bnf040302_maoi",
                        "040303"="bnf040303_ssri",
                        "040304"="bnf040304_oth"),
	   prescribed=1) %>%
distinct(eid, bnf, prescribed) %>%
pivot_wider(id_cols=eid, names_from=bnf, values_from=prescribed, values_fill=0) %>%
right_join(eids_with_records, by='eid') %>%
mutate(across(starts_with('bnf'), ~if_else(is.na(.), true=0, false=.))) %>%
mutate(bnf0403=if_any(starts_with('bnf'), ~ .x == 1)) %>%
mutate(bnf0403=bnf0403+0) %>%
mutate(FID=eid, IID=eid) %>%
select(FID, IID, starts_with('bnf')) |>
arrange(FID)


dir.create(here::here('data/phenotypes/health/health_outcomes/gp_scripts'))

write_tsv(gp_prescribed_ad, here::here('data/phenotypes/health/health_outcomes/gp_scripts/bnf_0403_antidepressants.tsv'))

# ATC codes

amp2_xml <- read_xml(here::here("data/phenotypes/health/TRUD/dmd/nhsbsa_dmd_7.2.0_20230717000001/f_amp2_3130723.xml"))

amp_xml <- xml_find_all(amp2_xml, ".//AMP")
amp_list <- as_list(amp_xml)

amp <- bind_rows(lapply(amp_list, unlist))

vpids_names <- amp |>
distinct(VPID, NM) |>
mutate(drug=str_to_upper(sapply(str_split(NM, pattern=" "), first))) |>
group_by(drug) |>
slice_head(n=1) |>
ungroup() |>
select(-NM)

bnf_xml <- read_xml(here::here("data/phenotypes/health/TRUD/dmdbonus/nhsbsa_dmdbonus_7.1.0_20230710000001/week282023-r2_3-BNF/f_bnf1_0060723.xml"))

bnf_list <- as_list(xml_find_first(bnf_xml, "VMPS"))

bnf_atc <- bind_rows(lapply(bnf_list, unlist)) |>
filter(!is.na(ATC), ATC != "n/a", !is.na(BNF)) 

bnf_atc_code <- bnf_atc |>
mutate(bnf_code = str_sub(BNF, 0, 7),
	  atc_code = str_sub(ATC, 0, 5)) |>
distinct(VPID, bnf_code, atc_code)

# ATC N06A
n06a <-
bnf_atc_code |>
filter(str_detect(atc_code, "^N06A")) |>
inner_join(vpids_names, by="VPID", multiple = "all")

gp_scripts_atc <-
gp_scripts_ad |>
mutate(drug=str_to_upper(sapply(str_split(drug_name, pattern=" "), first))) |>
left_join(n06a, by = "drug") |>
mutate(atc_code =
	case_when(drug %in% c("MOTIPRESS", "TRIPTAFEN") ~ "N06AA",
			 drug %in% c("BOLVIDON", "REBOXETINE", "VILOXAZINE", "VIVALAN", "VORTIOXETINE") ~ "N06AX",
			 TRUE ~ atc_code)
 ) |>
 filter(!is.na(atc_code)) |>
 distinct(eid, atc_code) |> 
 arrange(atc_code) |>
 mutate(value = 1) |>
 pivot_wider(names_from = "atc_code", values_from = "value", values_fill = 0)

atc_n0a6 <- 
eids_with_records |>
left_join(gp_scripts_atc, by = "eid") |>
mutate(across(starts_with("N06A"), ~ coalesce(.x, 0))) |>
mutate(N06A = if_any(starts_with("N06A"), ~ .x == 1)) |>
mutate(N06A = as.numeric(N06A)) |>
transmute(FID=eid, IID=eid, N06AA, N06AB, N06AX, N06A)

write_tsv(atc_n0a6, here::here('data/phenotypes/health/health_outcomes/gp_scripts/atc_N0A6_antidepressants.tsv'))