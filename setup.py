from setuptools import setup, find_packages
from pathlib import Path

here = Path(__file__).parent

about = {}
with open(here / "src" / "SimPG" / "__version__.py", encoding="utf-8") as f:
    exec(f.read(), about)


setup(
    name=about["__name__"],
    version=about["__version__"],
    author=about["__author__"],
    author_email=about["__author_email__"],
    description=about["__description__"],
    long_description=(here / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    url=about["__url__"],
    license=about["__license__"],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=["networkx==3.5"],
    python_requires=">=3.9",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    entry_points={
        "console_scripts": [
            "SimPG = SimPG.__main__:cli",
        ],
    },
)
