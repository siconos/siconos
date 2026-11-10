
This directory contains a set of options files that can be used to configure Siconos.


Usage:

- Copy or use one of the file in this directory, say somefile.cmake
- Edit its content to set options to suit your needs
- Then run


```bash

cmake -DUSER_OPTIONS_FILE=path_to/somefile.cmake ...
```




[default.cmake](./default.cmake) is the default when USER_OPTIONS_FILE is not set.

Each component can also be enabled or disabled directly from the command line
with a `WITH_<component>` option, for example:

```bash
cmake -S path_to_sources -B build-siconos -DWITH_numerics=ON -DWITH_kernel=ON
```

The available components are: `externals`, `numerics`, `kernel`, `control`,
`mechanics`, `mechanisms`, `io`, `vkernel`.  Dependencies are resolved
automatically, so enabling `mechanics` also enables `kernel`, `numerics` and
`externals`.

The legacy `COMPONENTS` variable is still supported: if it is set, it is used
as the source of truth and `WITH_<component>` options are ignored.


:warning: all files in the present directory are used in the Gitlab continuous integration process.

Please modify them with care (or not at all!)



Some extra options (designed for developers or advanced users) are available in the file [advanced_options.cmake](../cmake/advanced_options.cmake)
and can be set with the command line

```
cmake -DOPTION_NAME=OPTION_VALUE ...
```
