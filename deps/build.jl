using BinDeps

# TODO: Allow using HSL
@BinDeps.setup

libmetis = library_dependency("libmetis")
udir = "metis-4.0.3"
metis_dirname = joinpath(BinDeps.depsdir(libmetis), "src", udir)
provides(Sources, URI("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz"), libmetis, unpacked_dir=udir)

# compile_hsl = "MA57_SOURCE" in keys(ENV)
# if compile_hsl
# libmetis = library_dependency("libmetis")
# udir = "metis-4.0.3"
# metis_dirname = joinpath(BinDeps.depsdir(libmetis), "src", udir)
# provides(Sources, URI("http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz"), libmetis, unpacked_dir=udir)
  
  # libma57 = library_dependency("libhsl_ma57")
  # maudir = "hsl_ma57-5.2.0"
  # ma57_dirname = joinpath(BinDeps.depsdir(libma57), "src", maudir)
# end

# Metis
# provides(SimpleBuild,
#          (@build_steps begin
#             # Download and untar
#             GetSources(libmetis)
#             @build_steps begin
#               ChangeDirectory(BinDeps.depsdir(libmetis))        # Possibly remove
#               CreateDirectory("src")
#               CreateDirectory("usr")
#               CreateDirectory("usr/lib")
#               `tar -zxf downloads/metis-4.0.3.tar.gz -C src/` # Remove this later
#             end
#             @build_steps begin
#               ChangeDirectory(algencan_dirname)
#               # Compile with Makefile and flags
#               `make CFLAGS="-O3 -fPIC" FFLAGS="-O3 -ffree-form -fPIC"`
#               # Produce a shared library on deps/usr/lib
#               `gcc -shared -o ../../usr/lib/libalgencan.so
#                     -Wl,--whole-archive lib/libalgencan.a
#                     -Wl,--no-whole-archive -lgfortran`
#             end
          # end), libmetis, os = :Linux)

# provides(SimpleBuild,
#          (@build_steps begin
#             GetSources(libmetis)
#             @build_steps begin
#               ChangeDirectory(BinDeps.depsdir(libmetis))        # Possibly remove
#               CreateDirectory("src")
#               CreateDirectory("usr")
#               CreateDirectory("usr/lib")
#               `tar -zxf downloads/metis-4.0.3.tar.gz -C src/` # Remove this later
#             end
#             @build_steps begin
#               ChangeDirectory(metis_dirname)
#               `make COPTIONS=-fPIC`
#             end
#           end), libmetis, os = :Linux)

# # HSL
# provides(SimpleBuild,
#          (@build_steps begin
#             @build_steps begin
#               CreateDirectory(ma57_dirname)
#               FileUnpacker(ENV["MA57_SOURCE"], ma57_dirname, "")
#             end

#           end), libma57, os = :Linux)

# Algencan
provides(SimpleBuild,
         (@build_steps begin
            # Download and untar
            GetSources(libmetis)
            @build_steps begin
              ChangeDirectory(BinDeps.depsdir(libmetis))        # Possibly remove
              CreateDirectory("src")
              CreateDirectory("usr")
              CreateDirectory("usr/lib")
              `tar -zxf downloads/metis-4.0.3.tar.gz -C src/` # Remove this later
            end
            @build_steps begin
              ChangeDirectory(metis_dirname)
              # Compile with Makefile and flags
              `make COPTIONS="-fPIC"`
              # `make CFLAGS="-O3 -fPIC" FFLAGS="-O3 -ffree-form -fPIC"`
              # # Produce a shared library on deps/usr/lib
              # `gcc -shared -o ../../usr/lib/libalgencan.so
              #       -Wl,--whole-archive lib/libalgencan.a
              #       -Wl,--no-whole-archive -lgfortran`
            end
          end), libmetis, os = :Linux)

# # Algencan
# provides(SimpleBuild,
#          (@build_steps begin
#             # Download and untar
#             GetSources(libalgencan)
#             @build_steps begin
#               ChangeDirectory(BinDeps.depsdir(libalgencan))        # Possibly remove
#               CreateDirectory("src")
#               CreateDirectory("usr")
#               CreateDirectory("usr/lib")
#               `tar -zxf downloads/algencan-3.1.1.tgz -C src/` # Remove this later
#             end
#             @build_steps begin
#               ChangeDirectory(algencan_dirname)
#               # Compile with Makefile and flags
#               `make CFLAGS="-O3 -fPIC" FFLAGS="-O3 -ffree-form -fPIC"`
#               # Produce a shared library on deps/usr/lib
#               `gcc -shared -o ../../usr/lib/libalgencan.so
#                     -Wl,--whole-archive lib/libalgencan.a
#                     -Wl,--no-whole-archive -lgfortran`
#             end
#           end), libalgencan, os = :Linux)

# # TODO: see if it is possible to merge most of this two recipes.
# # This is mostly a dirty trick to get it compiling in OS X
# # as it does not accept --whole-archive in ld
# provides(SimpleBuild,
#          (@build_steps begin
#             # Download and untar
#             GetSources(libalgencan)
#             @build_steps begin
#               ChangeDirectory(BinDeps.depsdir(libalgencan))        # Possibly remove
#               CreateDirectory("src")
#               CreateDirectory("usr")
#               CreateDirectory("usr/lib")
#               `tar -zxf downloads/algencan-3.1.1.tgz -C src/` # Remove this later
#             end
#             @build_steps begin
#               ChangeDirectory(algencan_dirname)
#               # Compile with Makefile and flags
#               `make CFLAGS="-O3 -fPIC" FFLAGS="-O3 -ffree-form -fPIC"`
#               # Produce a shared library on deps/usr/lib
#               `gfortran -shared -o ../../usr/lib/libalgencan.dylib
#                     -Wl,-all_load lib/libalgencan.a
#                     -Wl,-noall_load -lgfortran`
#             end
#           end), libalgencan, os = :Darwin)

# @BinDeps.install Dict(:libalgencan => :libalgencan)
@BinDeps.install Dict(:libmetis => :libmetis)
