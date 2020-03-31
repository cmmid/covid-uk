library(readxl)
library(data.table)
library(socialmixr)
library(wpp2019)
library(countrycode)
library(ggplot2)

# POPULATIONS
data(pop)
data(UNlocations)
popF = data.table(popF)
popM = data.table(popM)
UNL = data.table(UNlocations)

codes = UNL[location_type <= 4, unique(country_code)]
populations = merge(popF[, c(1,2,3,18)], popM[, c(1,2,3,18)], by = c("country_code", "name", "age"))
names(populations)[4:5] = c("f", "m")
populations = populations[country_code %in% codes];
populations[, age := factor(age, levels = limits_to_agegroups(seq(0, 100, by = 5)))]
populations = merge(populations, UNL[, .(country_code, location_type)], by = "country_code")
populations = populations[order(location_type, name, age)]

# add China regional demographics. note populations expected in thousands, not in units
china_locs = c(
    Anhui = "^34",
    Guangdong = "^44",
    Guangxi = "^45",
    Hubei = "^42",
    Hunan = "^43",
    Jiangsu = "^32",
    Jiangxi = "^36",
    Jilin = "^22",
    Shaanxi = "^61",
    Shandong = "^37",
    Sichuan = "^51",
    Tianjin = "^12",
    Zhejiang = "^33",
    Wuhan = "^420100$", # Wuhan prefecture
    Beijing = "^110100$",
    Shanghai = "^310100$"
);

demog_china = data.table(readRDS("~/Dropbox/nCov-2019/data_sources/demographic_data/agedist_prf.rds")$Y2016);

for (i in 1:length(china_locs)) {
    adm_name = names(china_locs)[i];
    adm_code = china_locs[i];
    demog = demog_china[CNTY_CODE %like% adm_code, .(tot = sum(tot)), by = .(age)]
    demog[, age := as.numeric(gsub("-|\\+", "", age))];
    demog = demog[order(age)];
    demog[, age := paste0(age, "-", age + 4)];
    demog[age == "85-89", age := "85+"];
    demog[, age := factor(age, levels = limits_to_agegroups(seq(0, 85, by = 5)))]
    demog = demog[order(age)]
    populations = rbind(populations,
        demog[, .(country_code = 1111, name = paste("China |", adm_name), age = age, f = tot/2000, m = tot/2000, location_type = 5)]
    );
}

# demog_B = demog_china[CNTY_CODE == "110100"]
# demog_B[, age := as.numeric(gsub("-|\\+", "", age))];
# demog_B = demog_B[order(age)];
# demog_B[, age := paste0(age, "-", age + 4)];
# demog_B[age == "85-89", age := "85+"];
# populations = rbind(populations, demog_B[, .(country_code = 1111, name = "Beijing", age = age, f = tot/2000, m = tot/2000, location_type = 5)]);
# 
# demog_S = demog_china[CNTY_CODE == "310100"]
# demog_S[, age := as.numeric(gsub("-|\\+", "", age))];
# demog_S = demog_S[order(age)];
# demog_S[, age := paste0(age, "-", age + 4)];
# demog_S[age == "85-89", age := "85+"];
# populations = rbind(populations, demog_S[, .(country_code = 1111, name = "Shanghai", age = age, f = tot/2000, m = tot/2000, location_type = 5)]);

# load UK regional data
popUK = data.table(read_excel("~/Dropbox/nCoV/ukmidyearestimates20182019ladcodes.xlsx", "MYE2-All", "A5:CQ435"));
popUK[, Name := paste0("UK | ", Name)];
popUK[, Code := NULL];
popUK[, Geography1 := NULL];
popUK[, `All ages` := NULL];
popUK = melt(popUK, id.vars = "Name", variable.name = "age", value.name = "pop", variable.factor = F);
popUK[, age := as.numeric(age)]
popUK[, age_group := age %/% 5]
popUK[, age := paste0(age_group * 5, "-", age_group * 5 + 4)]
popUK[age == "90-94", age := "90+"]
popUK[, age := factor(age, levels = limits_to_agegroups(seq(0, 90, by = 5)))]
popUK = popUK[, .(f = sum(pop) / 2000, m = sum(pop) / 2000, location_type = 5, country_code = 2222), keyby = .(name = Name, age = age)]
populations = rbind(populations, popUK, fill = T)

# Bulawayo province
# NOTE CONTROVERSY OVER POPULATION SIZE
populations = rbind(populations, data.table(country_code = 3333, name = "Zimbabwe | Bulawayo",
    age = factor(limits_to_agegroups(seq(0, 85, by = 5)), levels = limits_to_agegroups(seq(0, 85, by = 5))),
    f = c(105116, 91354, 90964, 71443, 53485, 49776, 41871, 34063, 24302, 16397, 18154, 13762, 12493, 9858, 7613, 5481, 3946, 3260) / 2000,
    m = c(105116, 91354, 90964, 71443, 53485, 49776, 41871, 34063, 24302, 16397, 18154, 13762, 12493, 9858, 7613, 5481, 3946, 3260) / 2000,
    location_type = 5))

# Italy populations
# Source: for regions e.g. https://www.citypopulation.de/en/italy/admin/16__puglia/,
# for Milan https://www.citypopulation.de/en/italy/lombardia/milano/015146__milano/
italy_pop_structure = fread("~/Dropbox/nCoV/italy-pop-structure.txt")
italy_pop_structure = melt(italy_pop_structure, id.vars = c("age_lower", "age_upper"), variable.name = "location", value.name = "pop")

regioni = unique(italy_pop_structure$location)
popIT = NULL;
for (reg in regioni) {
    milan_target = italy_pop_structure[location == reg, pop/1000]
    milanpop = populations[name == "Italy"]
    milanpop[, name := paste("Italy |", reg)]
    milanpop[, location_type := 5]
    for (i in seq(1, 15, by = 2)) {
        targ = milan_target[i %/% 2 + 1];
        act = milanpop[i:(i+1), sum(f+m)]
        milanpop[i:(i+1), f := f * targ/act]
        milanpop[i:(i+1), m := m * targ/act]
    }

    targ = milan_target[9];
    act = milanpop[17:.N, sum(f+m)]
    milanpop[17:.N, f := f * targ/act]
    milanpop[17:.N, m := m * targ/act]

    # do a bit of smoothing
    for (i in seq(1, 15, by = 2)) {
        delta = (2 * milanpop[i + 1, .(f, m)] - milanpop[i, .(f, m)] - milanpop[i + 2, .(f, m)]) / 3
        milanpop[i, f := f + delta$f]
        milanpop[i, m := m + delta$m]
        milanpop[i + 1, f := f - delta$f]
        milanpop[i + 1, m := m - delta$m]
    }
    milanpop[, country_code := 5555]

    popIT = rbind(popIT, milanpop)
}

populations = rbind(populations, popIT)

# Certain major conurbations in LMIC
# Source: CIA world factbook, includes major urban areas in surrounding.
popLMIC = rbind(
    populations[name == "India",        .(country_code = 7777, name = "India | Delhi",               age, f = 30291 * f/sum(f+m), m = 30291 * m/sum(f+m), location_type)],
    populations[name == "Ethiopia",     .(country_code = 7777, name = "Ethiopia | Addis Ababa",      age, f =  4794 * f/sum(f+m), m =  4794 * m/sum(f+m), location_type)],
    populations[name == "Kenya",        .(country_code = 7777, name = "Kenya | Nairobi",             age, f =  4735 * f/sum(f+m), m =  4735 * m/sum(f+m), location_type)],
    populations[name == "South Africa", .(country_code = 7777, name = "South Africa | Johannesburg", age, f =  9677 * f/sum(f+m), m =  9677 * m/sum(f+m), location_type)],
    populations[name == "Nigeria",      .(country_code = 7777, name = "Nigeria | Lagos",             age, f = 14368 * f/sum(f+m), m = 14368 * m/sum(f+m), location_type)],
    populations[name == "Pakistan",     .(country_code = 7777, name = "Pakistan | Karachi",          age, f = 16094 * f/sum(f+m), m = 16094 * m/sum(f+m), location_type)],
    populations[name == "Bangladesh",   .(country_code = 7777, name = "Bangladesh | Dhaka",          age, f = 21006 * f/sum(f+m), m = 21006 * m/sum(f+m), location_type)]
);
populations = rbind(populations, popLMIC)

saveRDS(populations, "~/Dropbox/nCoV/wpp2019_pop2020.rds")




# CONTACT MATRICES

# Leung et al (HK)
hk_survey = get_survey("https://doi.org/10.5281/zenodo.1165561")

hk_pop = populations[name %like% "Hong Kong", f + m, by = age]$V1 * 1000
hk_survey_pop = data.frame(lower.age.limit = seq(0, 85, by = 5),
    population = c(hk_pop[1:17], sum(hk_pop[18:21])))

contact_matrix(hk_survey, age.limits = seq(0, 85, by = 5), survey.pop = hk_survey_pop, symmetric = T)$matrix

hk_home = contact_matrix(hk_survey, age.limits = seq(0, 85, by = 5), survey.pop = hk_survey_pop, filter = list(cnt_home = 1), symmetric = T)$matrix
hk_work = contact_matrix(hk_survey, age.limits = seq(0, 85, by = 5), survey.pop = hk_survey_pop, filter = list(cnt_work = 1), symmetric = T)$matrix
hk_scho = contact_matrix(hk_survey, age.limits = seq(0, 85, by = 5), survey.pop = hk_survey_pop, filter = list(cnt_school = 1), symmetric = T)$matrix
hk_othe = contact_matrix(hk_survey, age.limits = seq(0, 85, by = 5), survey.pop = hk_survey_pop, filter = list(cnt_home = 0, cnt_work = 0, cnt_school = 0), symmetric = T)$matrix

rownames(hk_home) = colnames(hk_home) = limits_to_agegroups(seq(0, 85, by = 5))
rownames(hk_work) = colnames(hk_work) = limits_to_agegroups(seq(0, 85, by = 5))
rownames(hk_scho) = colnames(hk_scho) = limits_to_agegroups(seq(0, 85, by = 5))
rownames(hk_othe) = colnames(hk_othe) = limits_to_agegroups(seq(0, 85, by = 5))

# Mossong et al (UK)
uk_survey = get_survey("https://doi.org/10.5281/zenodo.1043437")

uk_pop = populations[name %like% "United Kingdom", f + m, by = age]$V1 * 1000
uk_survey_pop = data.frame(lower.age.limit = seq(0, 75, by = 5),
    population = c(uk_pop[1:15], sum(uk_pop[16:21])))

uk_home = contact_matrix(uk_survey, countries = "United Kingdom", age.limits = seq(0, 75, by = 5), survey.pop = uk_survey_pop, filter = list(cnt_home = 1), symmetric = T)$matrix
uk_work = contact_matrix(uk_survey, countries = "United Kingdom", age.limits = seq(0, 75, by = 5), survey.pop = uk_survey_pop, filter = list(cnt_work = 1), symmetric = T)$matrix
uk_scho = contact_matrix(uk_survey, countries = "United Kingdom", age.limits = seq(0, 75, by = 5), survey.pop = uk_survey_pop, filter = list(cnt_school = 1), symmetric = T)$matrix
uk_othe = contact_matrix(uk_survey, countries = "United Kingdom", age.limits = seq(0, 75, by = 5), survey.pop = uk_survey_pop, filter = list(cnt_home = 0, cnt_work = 0, cnt_school = 0), symmetric = T)$matrix

rownames(uk_home) = colnames(uk_home) = limits_to_agegroups(seq(0, 75, by = 5))
rownames(uk_work) = colnames(uk_work) = limits_to_agegroups(seq(0, 75, by = 5))
rownames(uk_scho) = colnames(uk_scho) = limits_to_agegroups(seq(0, 75, by = 5))
rownames(uk_othe) = colnames(uk_othe) = limits_to_agegroups(seq(0, 75, by = 5))

# Mossong et al (Italy)
uk_survey = get_survey("https://doi.org/10.5281/zenodo.1043437")

it_pop = populations[name %like% "Italy", f + m, by = age]$V1 * 1000
it_survey_pop = data.frame(lower.age.limit = seq(0, 75, by = 5),
    population = c(uk_pop[1:15], sum(uk_pop[16:21])))

it_home = contact_matrix(uk_survey, countries = "Italy", age.limits = seq(0, 75, by = 5), survey.pop = uk_survey_pop, filter = list(cnt_home = 1), symmetric = T)$matrix
it_work = contact_matrix(uk_survey, countries = "Italy", age.limits = seq(0, 75, by = 5), survey.pop = uk_survey_pop, filter = list(cnt_work = 1), symmetric = T)$matrix
it_scho = contact_matrix(uk_survey, countries = "Italy", age.limits = seq(0, 75, by = 5), survey.pop = uk_survey_pop, filter = list(cnt_school = 1), symmetric = T)$matrix
it_othe = contact_matrix(uk_survey, countries = "Italy", age.limits = seq(0, 75, by = 5), survey.pop = uk_survey_pop, filter = list(cnt_home = 0, cnt_work = 0, cnt_school = 0), symmetric = T)$matrix

rownames(it_home) = colnames(it_home) = limits_to_agegroups(seq(0, 75, by = 5))
rownames(it_work) = colnames(it_work) = limits_to_agegroups(seq(0, 75, by = 5))
rownames(it_scho) = colnames(it_scho) = limits_to_agegroups(seq(0, 75, by = 5))
rownames(it_othe) = colnames(it_othe) = limits_to_agegroups(seq(0, 75, by = 5))


# Melegaro et al (Zimbabwe)
z_survey = get_survey("https://doi.org/10.5281/zenodo.1127693")
z_survey$contacts[, cnt_age_exact := as.character(cnt_age_exact)]
z_survey$contacts[cnt_age_exact %like% "[a-z]", cnt_age_exact := 0.5] # convert ages given in days, weeks, months to 0.5
z_survey$contacts[, cnt_age_exact := as.numeric(cnt_age_exact)]

# sample approxmate ages...
#z_survey$contacts[is.na(cnt_age_exact), cnt_age_exact := round(runif(1, cnt_age_est_min, cnt_age_est_max)), by = cont_id]

z_pop = populations[name %like% "Zimbabwe", f + m, by = age]$V1 * 1000
z_survey_pop = data.frame(lower.age.limit = seq(0, 85, by = 5),
    population = c(z_pop[1:17], sum(z_pop[18:21])))

z_home = contact_matrix(z_survey, age.limits = seq(0, 85, by = 5), survey.pop = z_survey_pop, filter = list(cnt_home = 1), symmetric = F, estimated.contact.age = "sample", missing.contact.age = "sample")$matrix
z_work = contact_matrix(z_survey, age.limits = seq(0, 85, by = 5), survey.pop = z_survey_pop, filter = list(cnt_work = 1), symmetric = F, estimated.contact.age = "sample", missing.contact.age = "sample")$matrix
z_scho = contact_matrix(z_survey, age.limits = seq(0, 85, by = 5), survey.pop = z_survey_pop, filter = list(cnt_school = 1), symmetric = F, estimated.contact.age = "sample", missing.contact.age = "sample")$matrix
z_othe = contact_matrix(z_survey, age.limits = seq(0, 85, by = 5), survey.pop = z_survey_pop, filter = list(cnt_home = 0, cnt_work = 0, cnt_school = 0), symmetric = F, estimated.contact.age = "sample", missing.contact.age = "sample")$matrix

rownames(z_home) = colnames(z_home) = limits_to_agegroups(seq(0, 85, by = 5))
rownames(z_work) = colnames(z_work) = limits_to_agegroups(seq(0, 85, by = 5))
rownames(z_scho) = colnames(z_scho) = limits_to_agegroups(seq(0, 85, by = 5))
rownames(z_othe) = colnames(z_othe) = limits_to_agegroups(seq(0, 85, by = 5))

matrices = list(
    "Hong Kong (Leung)" = 
        list(
            home = hk_home,
            work = hk_work,
            school = hk_scho,
            other = hk_othe
        ),
    "United Kingdom (Mossong)" =
        list(
            home = uk_home,
            work = uk_work,
            school = uk_scho,
            other = uk_othe
        ),
    "Italy (Mossong)" =
        list(
            home = it_home,
            work = it_work,
            school = it_scho,
            other = it_othe
        ),
    "Zimbabwe (Melegaro)" =
        list(
            home = z_home,
            work = z_work,
            school = z_scho,
            other = z_othe
        )
)

# Prem et al
read_prem_cm = function(n, country, loc = "all_locations")
{
    filename = paste0("~/Dropbox/nCoV/Readings/contact_matrices_152_countries/MUestimates_", loc, "_", n, ".xlsx");
    cm = as.matrix(read_excel(filename, country, col_names = ifelse(n == 1, T, F)));
    rownames(cm) = colnames(cm) = limits_to_agegroups(seq(0, 75, by = 5));
    return (cm)
}

# Get contact matrices
sheets1 = excel_sheets("~/Dropbox/nCoV/Readings/contact_matrices_152_countries/MUestimates_all_locations_1.xlsx");
sheets2 = excel_sheets("~/Dropbox/nCoV/Readings/contact_matrices_152_countries/MUestimates_all_locations_2.xlsx");


for (s in sheets1) {
    mh = read_prem_cm(1, s, "home");
    mw = read_prem_cm(1, s, "work");
    ms = read_prem_cm(1, s, "school");
    mo = read_prem_cm(1, s, "other_locations");
    
    matrices[[s]] = list(home = mh, work = mw, school = ms, other = mo);
}

for (s in sheets2) {
    mh = read_prem_cm(2, s, "home");
    mw = read_prem_cm(2, s, "work");
    ms = read_prem_cm(2, s, "school");
    mo = read_prem_cm(2, s, "other_locations");
    
    matrices[[s]] = list(home = mh, work = mw, school = ms, other = mo);
}




# getting group contacts from Zhang
part = merge(fread("~/Dropbox/nCoV/Contact Data/2019_Zhang_China_participant_common.csv"), fread("~/Dropbox/nCoV/Contact Data/2019_Zhang_China_participant_extra.csv"), by = "part_id");
cont = merge(fread("~/Dropbox/nCoV/Contact Data/2019_Zhang_China_contact_common.csv"), part[, .(part_id, part_age)], by = "part_id")
NROW = 18;

# group = function(age)
# {
#     pmin(age %/% 5 + 1, NROW)
# }
# 
# get.individual.contacts = function(part, cont, which)
# {
#     matk = matrix(0, nrow = NROW, ncol = NROW);
#     matu = matrix(0, nrow = NROW, ncol = NROW);
#     c_known = cont[!is.na(cnt_age_exact)];
#     
#     if (which == "h") {
#         c_known = c_known[cnt_home == 1];
#     } else if (which == "w") {
#         c_known = c_known[cnt_work == 1];
#     } else if (which == "s") {
#         c_known = c_known[cnt_school == 1];
#     } else {
#         c_known = c_known[cnt_home == 0 & cnt_work == 0 & cnt_school == 0];
#     }
#     
#     c_known[is.na(cnt_age_est_min), cnt_age_est_min := cnt_age_exact];
#     c_known[is.na(cnt_age_est_max), cnt_age_est_max := cnt_age_exact];
#     c_known[is.na(cnt_age_est_min), cnt_age_est_min := 0];
#     c_known[is.na(cnt_age_est_max), cnt_age_est_max := NROW*5];
#     c_unknown = cont[is.na(cnt_age_exact)];
#     
#     part_age_counts = part[, .N, keyby = group(part_age)]$N
#     
#     # do known contacts
#     for (r in 1:nrow(c_known))
#     {
#         pg = group(c_known[r, part_age]);
#         cg = group(c_known[r, runif(1, min(cnt_age_est_min, cnt_age_est_max), max(cnt_age_est_min, cnt_age_est_max))]);
#         matk[pg, cg] = matk[pg, cg] + 1 / part_age_counts[pg];
#     }
#     
#     # do unknown contacts
#     for (r in 1:nrow(c_unknown))
#     {
#         pg = group(c_unknown[r, part_age]);
#         if (sum(matk[pg,]) > 0) {
#             matu[pg,] = matu[pg,] + rmultinom(1, 1, matk[pg,]) / part_age_counts[pg];
#         } else {
#             lower = max(0, pg-1);
#             upper = min(NROW, pg+1);
#             matu[pg,] = matu[pg,] + rmultinom(1, 1, colMeans(matk[lower:upper,])) / part_age_counts[pg];
#         }
#     }
#     
#     matu + matk
# }
# 
# get.group.contacts = function(part, mat)
# {
#     grc = part[group_yn == 1]
#     mat = matrix(0, nrow = NROW, ncol = NROW)
#     
#     part_age_counts = part[, .N, keyby = group(part_age)]$N
#     
#     for (r in 1:nrow(grc))
#     {
#         pg = group(grc[r, part_age]);
#         n_cont_groups = grc[r, group_age0 + group_age6 + group_age12 + group_age19 + group_age60]
#         if (n_cont_groups == 0) {
#             possible_ages = 0:(NROW*5-1);
#         } else {
#             possible_ages = grc[r, c(rep(group_age0, 6), rep(group_age6, 6), rep(group_age12, 7), rep(group_age19, 41), rep(group_age60, (NROW*5)-60))]
#         }
# 
#         n_each = c(rmultinom(1, grc[r, group_n], possible_ages))
#         cgroups = group(rep(0:(NROW*5-1), n_each))
#         for (cg in cgroups) {
#             mat[pg, cg] = mat[pg, cg] + 1 / part_age_counts[pg];
#         }
#     }
#     mat
# }
# 
# nsamp = 50;
# 
# mat_h = matrix(0, nrow = NROW, ncol = NROW)
# mat_w = matrix(0, nrow = NROW, ncol = NROW)
# mat_s = matrix(0, nrow = NROW, ncol = NROW)
# mat_o = matrix(0, nrow = NROW, ncol = NROW)
# mat_g = matrix(0, nrow = NROW, ncol = NROW)
# 
# for (i in 1:nsamp) {
#     print(i)
#     mat_h = mat_h + get.individual.contacts(part, cont, "h") / nsamp
#     mat_w = mat_w + get.individual.contacts(part, cont, "w") / nsamp
#     mat_s = mat_s + get.individual.contacts(part, cont, "s") / nsamp
#     mat_o = mat_o + get.individual.contacts(part, cont, "o") / nsamp
#     mat_g = mat_g + get.group.contacts(part, mat) / nsamp
# }
# 
# image(mat_h+mat_w+mat_o+0.5*mat_g)
# 
# saveRDS(mat_h, "~/Dropbox/nCoV/shanghai_raw_h.rds")
# saveRDS(mat_w, "~/Dropbox/nCoV/shanghai_raw_w.rds")
# saveRDS(mat_s, "~/Dropbox/nCoV/shanghai_raw_s.rds")
# saveRDS(mat_o, "~/Dropbox/nCoV/shanghai_raw_o.rds")
# saveRDS(mat_g, "~/Dropbox/nCoV/shanghai_raw_g.rds")

mat_h = readRDS("~/Dropbox/nCoV/shanghai_raw_h.rds")
mat_w = readRDS("~/Dropbox/nCoV/shanghai_raw_w.rds")
mat_s = readRDS("~/Dropbox/nCoV/shanghai_raw_s.rds")
mat_o = readRDS("~/Dropbox/nCoV/shanghai_raw_o.rds")
mat_g = readRDS("~/Dropbox/nCoV/shanghai_raw_g.rds")

symmetrize = function(mat, pop)
{
    norm.mat = diag(pop) %*% mat;
    return (0.5 * diag(1/pop) %*% (norm.mat + t(norm.mat)));
}

add_names = function(mat)
{
    dimnames(mat) = list(
        c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+"),
        c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
    )
    # dimnames(mat) = list(
    #     c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+"),
    #     c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+")
    # )
    mat
}

# plUK = sp_contact_matrices("United Kingdom (Mossong)", matrices[["United Kingdom (Mossong)"]], legpos = "right")

# # ASSUMPTION 1 - all group contacts involving 0-19 year olds are school
# mat_g_school = mat_g
# mat_g_nonschool = mat_g
# mat_g_school[5:18, 5:18] = 0
# mat_g_nonschool[, 1:4] = 0
# mat_g_nonschool[1:4, ] = 0
# 
# wuhan_pop = populations[name == "Wuhan", f + m]
# 
# wuhan_new_matrices = list(
#     home   = add_names(symmetrize(mat_h, wuhan_pop)),
#     work   = add_names(symmetrize(mat_w, wuhan_pop)),
#     school = add_names(symmetrize(mat_s + mat_g_school, wuhan_pop)),
#     other  = add_names(symmetrize(mat_o + mat_g_nonschool, wuhan_pop))
# )
# 
# # check
# with(matrices[["United Kingdom (Mossong)"]], rowSums(school)/rowSums(work + home + other))
# with(wuhan_new_matrices, rowSums(school)/rowSums(work + home + other))
# 
# plWu1 = sp_contact_matrices("Wuhan (group = full strength; school is all with 0-19)", wuhan_new_matrices, legpos = "right")

# ASSUMPTION 2 - only contacts between 0-19 year olds are school, and group contacts are worth 0.5 as much as individual contacts
mat_g_school = mat_g
mat_g_nonschool = mat_g
mat_g_school[1:4, 5:NROW] = 0
mat_g_school[5:NROW, 1:NROW] = 0
mat_g_nonschool[1:4, 1:4] = 0

china_pops = populations[name %like% "^China \\| ", unique(name)]
for (reg_name in china_pops) {
    reg_pop = populations[name == reg_name, f + m];
    reg_pop = c(reg_pop[1:NROW-1], sum(reg_pop[NROW:length(reg_pop)]));

    reg_matrices = list(
        home   = add_names(symmetrize(mat_h, reg_pop)),
        work   = add_names(symmetrize(mat_w, reg_pop)),
        school = add_names(symmetrize(mat_s + 0.5*mat_g_school, reg_pop)),
        other  = add_names(symmetrize(mat_o + 0.5*mat_g_nonschool, reg_pop))
    )

    matrices[[reg_name]] = reg_matrices
}

# All China
C_pop = populations[name == "China", f + m]
C_pop = c(C_pop[1:(NROW-1)], sum(C_pop[NROW:length(C_pop)]))

C_new_matrices = list(
    home   = add_names(symmetrize(mat_h, C_pop)),
    work   = add_names(symmetrize(mat_w, C_pop)),
    school = add_names(symmetrize(mat_s + 0.5*mat_g_school, C_pop)),
    other  = add_names(symmetrize(mat_o + 0.5*mat_g_nonschool, C_pop))
)

matrices[["China | China"]] = C_new_matrices

# UK regional contact matrices
UKregions = populations[name %like% "^UK \\| ", unique(name)]

for (reg in UKregions) {
    print(reg)
    
    ruk_pop = populations[name == reg, f + m, by = age]$V1 * 1000
    ruk_survey_pop = data.frame(lower.age.limit = seq(0, 75, by = 5),
        population = c(ruk_pop[1:15], sum(ruk_pop[16:19])))

    ruk_home = contact_matrix(uk_survey, countries = "United Kingdom", age.limits = seq(0, 75, by = 5), 
        survey.pop = ruk_survey_pop, filter = list(cnt_home = 1), symmetric = T)$matrix
    ruk_work = contact_matrix(uk_survey, countries = "United Kingdom", age.limits = seq(0, 75, by = 5), 
        survey.pop = ruk_survey_pop, filter = list(cnt_work = 1), symmetric = T)$matrix
    ruk_scho = contact_matrix(uk_survey, countries = "United Kingdom", age.limits = seq(0, 75, by = 5), 
        survey.pop = ruk_survey_pop, filter = list(cnt_school = 1), symmetric = T)$matrix
    ruk_othe = contact_matrix(uk_survey, countries = "United Kingdom", age.limits = seq(0, 75, by = 5), 
        survey.pop = ruk_survey_pop, filter = list(cnt_home = 0, cnt_work = 0, cnt_school = 0), symmetric = T)$matrix

    rownames(ruk_home) = colnames(ruk_home) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(ruk_work) = colnames(ruk_work) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(ruk_scho) = colnames(ruk_scho) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(ruk_othe) = colnames(ruk_othe) = limits_to_agegroups(seq(0, 75, by = 5))
    
    matrices[[reg]] = list(
        home = ruk_home,
        work = ruk_work,
        school = ruk_scho,
        other = ruk_othe
    );
}

# Italy regional contact matrices
regioni = populations[name %like% "^Italy \\| ", unique(name)]

for (reg in regioni) {
    print(reg)
    
    rit_pop = populations[name == reg, f + m, by = age]$V1 * 1000
    rit_survey_pop = data.frame(lower.age.limit = seq(0, 75, by = 5),
        population = c(rit_pop[1:15], sum(rit_pop[16:19])))

    rit_home = contact_matrix(uk_survey, countries = "Italy", age.limits = seq(0, 75, by = 5), 
        survey.pop = rit_survey_pop, filter = list(cnt_home = 1), symmetric = T)$matrix
    rit_work = contact_matrix(uk_survey, countries = "Italy", age.limits = seq(0, 75, by = 5), 
        survey.pop = rit_survey_pop, filter = list(cnt_work = 1), symmetric = T)$matrix
    rit_scho = contact_matrix(uk_survey, countries = "Italy", age.limits = seq(0, 75, by = 5), 
        survey.pop = rit_survey_pop, filter = list(cnt_school = 1), symmetric = T)$matrix
    rit_othe = contact_matrix(uk_survey, countries = "Italy", age.limits = seq(0, 75, by = 5), 
        survey.pop = rit_survey_pop, filter = list(cnt_home = 0, cnt_work = 0, cnt_school = 0), symmetric = T)$matrix

    rownames(rit_home) = colnames(rit_home) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(rit_work) = colnames(rit_work) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(rit_scho) = colnames(rit_scho) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(rit_othe) = colnames(rit_othe) = limits_to_agegroups(seq(0, 75, by = 5))
    
    matrices[[reg]] = list(
        home = rit_home,
        work = rit_work,
        school = rit_scho,
        other = rit_othe
    );
}

# Zimbabwe regional contact matrices
zimreg = populations[name %like% "^Zimbabwe \\| ", unique(name)]

msum = function(mat)
{
    m = mat$matrices[[1]]$matrix * 1/length(mat$matrices)
    for (i in 2:length(mat$matrices)) {
        m = m + mat$matrices[[i]]$matrix * 1/length(mat$matrices)
    }
    m
}

for (reg in zimreg) {
    print(reg)
    
    rzm_pop = populations[name == reg, f + m, by = age]$V1 * 1000
    rzm_survey_pop = data.frame(lower.age.limit = seq(0, 75, by = 5),
        population = c(rzm_pop[1:15], sum(rzm_pop[16:18])))

    rzm_home = msum(contact_matrix(z_survey, age.limits = seq(0, 75, by = 5), 
        survey.pop = rzm_survey_pop, filter = list(cnt_home = 1), symmetric = T, estimated.contact.age = "sample", n = 10))
    rzm_work = msum(contact_matrix(z_survey, age.limits = seq(0, 75, by = 5), 
        survey.pop = rzm_survey_pop, filter = list(cnt_work = 1), symmetric = T, estimated.contact.age = "sample", n = 10))
    rzm_scho = msum(contact_matrix(z_survey, age.limits = seq(0, 75, by = 5), 
        survey.pop = rzm_survey_pop, filter = list(cnt_school = 1), symmetric = T, estimated.contact.age = "sample", n = 10))
    rzm_othe = msum(contact_matrix(z_survey, age.limits = seq(0, 75, by = 5), 
        survey.pop = rzm_survey_pop, filter = list(cnt_home = 0, cnt_work = 0, cnt_school = 0), symmetric = T, estimated.contact.age = "sample", n = 10))

    rownames(rzm_home) = colnames(rzm_home) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(rzm_work) = colnames(rzm_work) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(rzm_scho) = colnames(rzm_scho) = limits_to_agegroups(seq(0, 75, by = 5))
    rownames(rzm_othe) = colnames(rzm_othe) = limits_to_agegroups(seq(0, 75, by = 5))
    
    matrices[[reg]] = list(
        home = rzm_home,
        work = rzm_work,
        school = rzm_scho,
        other = rzm_othe
    );
}

saveRDS(matrices, "~/Dropbox/nCoV/all_matrices.rds");
