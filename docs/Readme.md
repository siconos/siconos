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

--> to create tex files


# Transfer to Gforge: 


scp -pr siconos_build_dir/build/docs/build/html/* login@scm.gforge.inria.fr:/home/groups/siconos/htdocs

-p option for first commit only => allow to keep same files permissions (rw|rw|r) after sending. 

Warning: if you generate latex or xml do not transfer it! 



# Tests: use gitlab to create and publish siconos doc

!! Only in update_doc branch

## Details:

* A mirror of siconos project on gricad-gitlab : https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-mirror
* CI is configured to build and publish documentation (see files .gitlab-ci.yml and CI/make_siconos_doc.sh)

	* gitlab-ci.yml : two jobs (make_doc and publish)
	* make_siconos_doc.sh : describe building process (apt + install python packages from requirements.txt + cmake and make doc)
	 Siconos conf is described in CI/siconos_docs.cmake


* Check CI status : https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-mirror/pipelines
* Check (public) website : https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-mirror


## Usage:

* Add remote repo:

```
git remote add gitlab_for_doc git@gricad-gitlab.univ-grenoble-alpes.fr:nonsmooth/siconos-mirror.git
```

* To publish doc:

```
git checkout update_doc
git push gitlab_for_doc update_doc
