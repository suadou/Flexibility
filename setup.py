
from distutils.core import setup


setup(name='flexprot',
      version='1.0',
      description='Assesment of protein flexibility',
      author='Silvia González López & Sergio Suárez Dou',
      author_email='silvia.gonzalez10@estudiant.upf.edu',
      packages=['my_package', 'my_package.subpackage'],
      install_requires=[
          'Bio',
          'scipy',
          'pandas',
          'matplotlib',
          'seaborn',
          'argparse'
          ])
