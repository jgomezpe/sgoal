import setuptools
from sgoal.version import Version


setuptools.setup(name='sgoal',
      version=Version('1.0.0').number,
      description='Stochastic Global Optimization Algorithms Python Package',
      long_description=open('README.md').read().strip(),
      author='Jonatan Gomez',
      author_email='jgomezpe@unal.edu.co',
      url='https://github.com/jgomezpe/sgoal',
      packages=setuptools.find_packages(),
      install_requires=[],
      license='MIT License',
      zip_safe=False,
      keywords='sgoal package',
      classifiers=['Packages', 'SGoal']
)