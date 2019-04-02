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

* A mirror of siconos project on gricad-gitlab : https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos
* CI is configured to build and publish documentation (see files .gitlab-ci.yml and CI/make_siconos_doc.sh)

	* gitlab-ci.yml : two jobs (make_doc and publish)
	* make_siconos_doc.sh : describe building process (apt + install python packages from requirements.txt + cmake and make doc)
	 Siconos conf is described in CI/siconos_docs.cmake

* Check CI status : https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-mirror/pipelines
* Check (public) website : https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-mirror


## Usage:

* Add remote repo:

```
git remote add sico-doc git@gricad-gitlab.univ-grenoble-alpes.fr:nonsmooth/siconos-mirror.git
```

* To publish doc:

```
# Update your repo with github master
git checkout master
git pull --rebase origin master git checkout master
# Push to gitlab (--> will automatically push to github master)
git push sico-doc master

