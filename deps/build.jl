using BinDeps

@BinDeps.setup

ma57_src = ENV["MA57_SOURCE"]
libhsl_ma57 = library_dependency("libhsl_ma57")
src_dir = joinpath(BinDeps.depsdir(libhsl_ma57), "src")
ma57_dir = joinpath(BinDeps.depsdir(libhsl_ma57), "src", "hsl_ma57-5.2.0")

# HSL
provides(SimpleBuild, 
  (@build_steps begin
      CreateDirectory(ma57_dir)
      `tar xvf $ma57_src --directory=$src_dir`
  end), libhsl_ma57, os = :Linux)

@BinDeps.install Dict(:libhsl_ma57 => :libhsl_ma57)
