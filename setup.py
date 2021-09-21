from setuptools import setup

setup(
   name='PropTools',
   version='0.1',
   description='PropTools',
   author='Tom Aveyard',
   author_email='thomasaveyard@gmail.com',
   packages=['PropTools'],
   install_requires=['ambiance', 'CoolProp', 'matplotlib', 'numpy', 'pandas', 'rocketcea', 'scipy', 'thermo']
)


