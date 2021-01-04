using Documenter, NLPModelsAlgencan

makedocs(
   sitename="NLPModelsAlgencan.jl",
   pages = [
      "Home" => "index.md",
      "Get started" => [
         "First steps" => "first_steps.md",
         "Optional parameters" => "parameters.md"
      ],
      "Examples of usage" => "examples.md"
   ]
)

deploydocs(
    repo="github.com/pjssilva/NLPModelsAlgencan.jl.git"
)
