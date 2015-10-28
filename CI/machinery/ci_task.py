class CiTask():

   def __init__(self,
                mode='Continuous',
                build_configuration='Release',
                distrib=None,
                ci_config=None,
                fast=False,
                pkgs=None):
      self._fast=fast
      self._distrib=distrib
      self._mode=mode
      self._build_configuration=build_configuration
      self._ci_config=ci_config
      self._pkgs=pkgs

   def templates(self):
      return ','.join(self._pkgs)

   def copy(self):
      def init(mode=self._mode, build_configuration=self._build_configuration,
               distrib=self._distrib,
               ci_config=self._ci_config, fast=self._fast, pkgs=self._pkgs,
               add_pkgs=None, remove_pkgs=None):
         if add_pkgs is not None:
            pkgs = self._pkgs + add_pkgs

         if remove_pkgs is not None:
            pkgs = filter(lambda p: p not in remove_pkgs, pkgs)

         return CiTask(mode, build_configuration, distrib, ci_config, fast,
                       pkgs)
      return init
