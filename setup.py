from setuptools import setup, Extension, find_packages


setup(
    name="ccmgen-scripts",
    version="1.0.0",
    description="Reproducing analysis of study for CCMgen",
    license="AGPLv3",
    author="Susann Vorberg, Stefan Seemayer, Johannes Soeding",
    author_email="Susann.Vorberg@gmail.com",
    url="https://github.com/soedinglab/ccmgen-scripts",
    packages=find_packages(),
    install_requires=['msgpack-python', 'numpy', 'plotly', 'colorlover']
)
