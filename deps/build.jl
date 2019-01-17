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
      @build_steps begin
        ChangeDirectory(ma57_dir)
        `patch -p1 <../../patches/patch_ma57.txt`
        `./configure --prefix=$ma57_dir CFLAGS=-fPIC FCFLAGS=-fPIC `
        `make`
        `make install`
      end
  end), libhsl_ma57, os = :Linux)

@BinDeps.install Dict(:libhsl_ma57 => :libhsl_ma57)
