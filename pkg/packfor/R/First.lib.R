.onAttach <- function(lib, pkg)  {
	packageStartupMessage("packfor: R Package for Forward Selection (Canoco Manual p.49)",appendLF = TRUE)
	packageStartupMessage("version",utils::packageDescription("packfor",field="Version"),appendLF = TRUE)
}
