import setuptools
version = "0.0.14"

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setuptools.setup(
    name="symlib",
    version=version,
    author="Philip Mansfield",
    author_email="mansfield.astro@gmail.com",
    description="A library for working with data from the Symphony and MW-est zoom-in suites.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=["numpy", "scipy", "matplotlib", "colossus"],
    url="https://github.com/phil-mansfield/symlib",
    keywords=["python"],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.6"
)
#print(setuptools.find_packages(where="src"))
