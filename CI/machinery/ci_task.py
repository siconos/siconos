class CiTask():

   def __init__(self,
                mode='Continuous',
                distrib='ubuntu:14.04',
                ci_config='default',
                cmake_args=[],
                pkgs=['build-base','gcc', 'gfortran', 'g++', 'atlas-lapack']):
      self._distrib=distrib
      self._mode=mode
      self._ci_config=ci_config
      self._cmake_args=cmake_args
      self._pkgs=pkgs

   def templates(self):
      return ','.join(self._pkgs)

