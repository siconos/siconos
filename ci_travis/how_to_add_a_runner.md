Process to add a new gitlab-ci runner for the nonsmooth group
-------------------------------------------------------------

On the remote host of the runner, as root

Change home dir of gitlab-runner user (since on inria's machine, /home is not writable)
```
useradd gitlab-runner -d /scratch/gitlab-runner
```

install docker and gitlab-runner and register

```
apt install docker.io
curl -L https://packages.gitlab.com/install/repositories/runner/gitlab-runner/script.deb.sh | sudo bash
apt install gitlab-runner
gitlab-runner register
```

Input values for register are given on https://gricad-gitlab.univ-grenoble-alpes.fr/groups/nonsmooth/-/settings/ci_cd.

Use 'docker' as executor and ruby:2.3 as default Docker image.

That's it.
The runner must appear in the list on https://gricad-gitlab.univ-grenoble-alpes.fr/groups/nonsmooth/-/settings/ci_cd
and is from now on available for all projects of the group.


