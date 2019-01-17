using BinDeps

@BinDeps.setup

ma57_src = ENV["MA57_SOURCE"]
libhsl_ma57 = library_dependency("libhsl_ma57")
ma57_dir = joinpath(BinDeps.depsdir(libhsl_ma57), "src", "hsl_ma57-5.2.0")
metis_dir = joinpath(BinDeps.depsdir(libhsl_ma57), "src", "metis-4.0.3")
provides(Sources, URI("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz"), libhsl_ma57, unpacked_dir="metis-4.0.3")

# src_dir = joinpath(BinDeps.depsdir(libhsl_ma57), "src")
src_dir = joinpath(BinDeps.depsdir(libhsl_ma57), "src")

# HSL
provides(SimpleBuild, 
  (@build_steps begin
      # Get Metis sources and unpack
      ChangeDirectory(inDeps.depsdir(libhsl_ma57))
      GetSources(libhsl_ma57)
      `tar xvf downloads/metis-4.0.3.tar.gz --directory=$src_dir`

      # Unpack HSL 
      CreateDirectory(ma57_dir)
      `tar xvf $ma57_src --directory=$src_dir`

      # Build HSL
      @build_steps begin
        ChangeDirectory(ma57_dir)
        `patch -p1 <../../patches/patch_ma57.txt`
        `./configure --prefix=$ma57_dir CFLAGS=-fPIC FCFLAGS=-fPIC `
        `make`
        `make install`
      end
      @build_steps begin
        ChangeDirectory(joinpath(ma57_dir, "lib"))
        `gcc --shared -o libhsl_ma57.so libhsl_ma57.a`
      end
  end), libhsl_ma57, os = :Linux
)

@BinDeps.install Dict(:libhsl_ma57 => :libhsl_ma57)
