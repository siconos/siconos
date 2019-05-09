The build-siconos-vs2013-64.bat file needs to be customized to run the build process

The post-checkout file is to be used as a git hook to deal with the symlinks in the git repository, since Windows does not support those
The deal-with-symlink.sh replace the dead symlink by the actual file
