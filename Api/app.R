library(plumber)
library(HMR)
library(dplyr)

Plumber <- plumber::plumb("plumber.R");
SPEC=Plumber$getApiSpec();

for(action in names(SPEC$paths)){
  writeLines(action)

  for(verb in names(SPEC$paths[[action]])){
    writeLines(verb)
    SPEC$paths[[action]][[verb]]$operationId=paste(action,verb,sep ="_");
  }
}

SPEC$paths$`/calculate`$post = within(SPEC$paths$`/calculate`$post, rm('parameters'))
#SPEC$paths$`/calculate`$post = within(SPEC$paths$`/calculate`$post, rm('requestBody'))
SPEC$paths$`/calculate`$post = within(SPEC$paths$`/calculate`$post, rm('responses'))

SPEC$paths$`/calculate`$post$requestBody$content$`application/json`$schema$type = "array"
SPEC$paths$`/calculate`$post$requestBody$content$`application/json`$schema$items$`$ref` = "#/components/schemas/ConcentrationMeasurement"

SPEC$paths$`/calculate`$post$responses$'200'$content$'application/json'$schema$type = "array"
SPEC$paths$`/calculate`$post$responses$'200'$content$'application/json'$schema$items$`$ref` = "#/components/schemas/CalculateResponse"
SPEC$paths$`/calculate`$post$responses$'200'$description = "OK"

SPEC$components$schemas$ConcentrationMeasurement$type = "object"
SPEC$components$schemas$ConcentrationMeasurement$required = list("SeriesIdentifier","GasConcentration","ChamberVolume","ChamberArea","DeploymentTime")
SPEC$components$schemas$ConcentrationMeasurement$properties$SeriesIdentifier$type = "string"
SPEC$components$schemas$ConcentrationMeasurement$properties$GasConcentration$type = "number"
SPEC$components$schemas$ConcentrationMeasurement$properties$GasConcentration$format = "double"
SPEC$components$schemas$ConcentrationMeasurement$properties$ChamberVolume$type = "number"
SPEC$components$schemas$ConcentrationMeasurement$properties$ChamberVolume$format = "double"
SPEC$components$schemas$ConcentrationMeasurement$properties$ChamberArea$type = "number"
SPEC$components$schemas$ConcentrationMeasurement$properties$ChamberArea$format = "double"
SPEC$components$schemas$ConcentrationMeasurement$properties$DeploymentTime$type = "number"
SPEC$components$schemas$ConcentrationMeasurement$properties$DeploymentTime$format = "double"

SPEC$components$schemas$CalculateResponse$type = "object"
SPEC$components$schemas$CalculateResponse$required = list("SeriesIdentifier","SuggestedMethod","HMRResult","LinearResult")

SPEC$components$schemas$CalculateResponse$properties$SeriesIdentifier$type = "string"

SPEC$components$schemas$CalculateResponse$properties$SuggestedMethod$type = "string"
SPEC$components$schemas$CalculateResponse$properties$SuggestedMethod$enum = list("Linear","HMR","NoFlux")
SPEC$components$schemas$CalculateResponse$properties$SuggestedMethod$`x-ms-enum`$name = "FluxMethod"
SPEC$components$schemas$CalculateResponse$properties$SuggestedMethod$`x-ms-enum`$modelAsString = FALSE
SPEC$components$schemas$CalculateResponse$properties$LinearResult$`$ref` = "#/components/schemas/LinearResultModel"
SPEC$components$schemas$CalculateResponse$properties$HMRResult$`$ref` = "#/components/schemas/HMRResultModel"


SPEC$components$schemas$CalculateResponse$properties$Warning$type = "string"

SPEC$components$schemas$ResultModelBase$type = "object"
SPEC$components$schemas$ResultModelBase$required = list("FluxValue","FluxStandardError","FluxPValue","FluxLower95","FluxUppoer95")
SPEC$components$schemas$ResultModelBase$properties$FluxValue$type = "number"
SPEC$components$schemas$ResultModelBase$properties$FluxValue$format = "double"
SPEC$components$schemas$ResultModelBase$properties$FluxStandardError$type = "number"
SPEC$components$schemas$ResultModelBase$properties$FluxStandardError$format = "double"
SPEC$components$schemas$ResultModelBase$properties$FluxPValue$type = "number"
SPEC$components$schemas$ResultModelBase$properties$FluxPValue$format = "double"
SPEC$components$schemas$ResultModelBase$properties$FluxLower95$type = "number"
SPEC$components$schemas$ResultModelBase$properties$FluxLower95$format = "double"
SPEC$components$schemas$ResultModelBase$properties$FluxUppoer95$type = "number"
SPEC$components$schemas$ResultModelBase$properties$FluxUppoer95$format = "double"

resultmodelbaseref = list(`$ref` = "#/components/schemas/ResultModelBase")

hmritem = list(type="object",required=list("KappaValue","HValue"),properties = list(HValue = list(type = "number",format="double"), KappaValue = list(type = "number",format="double")))
hmritem = list(resultmodelbaseref,hmritem)
SPEC$components$schemas$HMRResultModel$allOf = hmritem

linearitem = list(type="object",required=list("InterceptValue"),properties = list(InterceptValue = list(type = "number",format="double")))
linearitem = list(resultmodelbaseref,linearitem)
SPEC$components$schemas$LinearResultModel$allOf = linearitem


Plumber$setApiSpec(SPEC);


Plumber$run(host = '0.0.0.0', port = 80,swagger = TRUE );