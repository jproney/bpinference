context("test conversion of real data to stan-friendly format")

test_that("Incorrect inputs trigger apropriate errors" , {
  lincs = read.csv("../../lincs-data/timeseries/lincs_ts.csv")
  lincs <- lincs %>% group_by(AssayWell, Plate) %>% mutate(PrevCellCount = dplyr::lag(CellCountAfterTreatment), DeltaT = Timepoint - dplyr::lag(Timepoint))
  lincs <- lincs %>% filter(!is.na(PrevCellCount))
  controls <- lincs %>% filter(is.na(Drug) & Timepoint < 40)
  
  ggplot(controls, aes(x = Timepoint, y = CellCountAfterTreatment)) + geom_point(aes(color = AssayWell))
  
  func_deps = c("c[1]","c[2]")
  ndep = 0
  nparam = 2
  mod = bp_model_simple_birth_death(func_deps, nparam, ndep)
  
  expect_error(create_stan_data(mod, controls$CellCountAfterTreatment, controls$PrevCellCount, controls$DeltaT[-1], simple_bd = T),"times, initial populations, and final populations must all have same number of rows!")
  expect_error(create_stan_data(mod, controls$CellCountAfterTreatment, controls$PrevCellCount[-1], controls$DeltaT, simple_bd = T),"times, initial populations, and final populations must all have same number of rows!")
  expect_error(create_stan_data(mod, controls$CellCountAfterTreatment[-1], controls$PrevCellCount, controls$DeltaT, simple_bd = T),"times, initial populations, and final populations must all have same number of rows!")
  expect_warning(create_stan_data(mod, controls$CellCountAfterTreatment, matrix(controls$PrevCellCount, ncol=1), controls$DeltaT, simple_bd = T),"final_pop is a vector, not a matrix. Converting to a one-column matrix")
  
})
