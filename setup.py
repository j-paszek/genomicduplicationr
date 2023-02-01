import setuptools

with open("README.md", "r", encoding="utf-8") as ld:
    long_description = ld.read()

setuptools.setup(
    name="gdscore",
    version="1.0",
    author="Jaroslaw Paszek",
    author_email="jpaszek@mimuw.edu.pl",
    license="Creative Commons Attribution 4.0 International",
    description="Tool for computing genomic duplications obtained by single gene duplication clustering. "
                "Available clustering methods - EC and ME. Available models - LCA, GMS, PG and FHS (unrestricted).",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/j-paszek/genomicduplicationr.git",
    packages=setuptools.find_packages(),
    # packages=['gdscore', 'tests'],
    entry_points={'console_scripts': ['gdscore=gdscore.gdscore:main']},
    python_requires='>=3.6',
    install_requires=['pytest==6.2.4']
)
