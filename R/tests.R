dat = readr::read_csv('R:/landscape_ecology/projects/winch/data/inclinometer_load_cell_merged.csv')
dat = dplyr::filter(dat, treeid=='JK018')
dat = dplyr::select(dat, F_kN, t_angle)

moment_kNm = dat$F_kN * 10
tilt_deg = dat$t_angle
diam_cm = 20
ht_m = 10

ss_df = stress_strain(moment_kNm, tilt_deg, diam_cm, ht_m)

getMOE(ss_df$stress, ss_df$strain)
