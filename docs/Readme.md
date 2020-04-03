# Web pages generation 

in siconos build directory:

```
cmake ... -DWITH_DOCUMENTATION=ON
make doc
```

Build html documentation in siconos/build_dir/docs/build/html

```
make latex
```

--> to create tex files (experimental !)


# To publish web site and documentation :

* just add [doc-build] in your commit message. This will automatically activate a job of continuous integration that will create and publish the html documentation.

* /!\ Push to gricad-gitlab project, not to github, the mirroring process will do the job /!\

* Check CI status: https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/pipelines

* Check (public) website : https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos

Note for developers:

the whole documentation concerning the documentation process is in the developers' manual,

https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/devel_guide/write_and_build_doc.html

