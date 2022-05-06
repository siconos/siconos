# Web pages generation 

in siconos build directory:

```
cmake ... -DWITH_DOCUMENTATION=ON
make doc
```

Build html documentation in siconos/build_dir/docs/build/html


# To publish web site and documentation :

* Visit CI status page of Siconos on gricad-gitlab:

	https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/pipelines
	
Pick a pipeline and "play" (manual) the generate-doc job.
	
* Check (public) website : https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos

Note for developers:

the whole documentation concerning the documentation process is in the developers' manual,

https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/devel_guide/write_and_build_doc.html

