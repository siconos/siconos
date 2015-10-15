class CiTask():

   def __init__(self,
                mode='Continuous',
                distrib='ubuntu:14.04',
                ci_config='default',
                pkgs=['siconos-default','gcc', 'gfortran', 'g++', 'atlas-lapack']):
      self._distrib=distrib
      self._mode=mode
      self._ci_config=ci_config
      self._pkgs=pkgs

   def templates(self):
      return ','.join(self._pkgs)
