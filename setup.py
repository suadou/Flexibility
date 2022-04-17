from setuptools import setup, find_packages

setup(name='flexprot',
      version='1.0',
      description='Assesment of protein flexibility',
      author='Silvia González López & Sergio Suárez Dou',
      author_email='silvia.gonzalez10@estudiant.upf.edu',
      packages=find_packages(),
      install_requires=[
          'biopython',
          'scipy',
          'pandas',
          'matplotlib',
          'seaborn',
          'argparse',
          'threaded',
          'configparser',
          'numpy',
          'regex',
          'python-math'
          ])
