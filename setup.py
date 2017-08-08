try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name = "index_bam_by_read_id",
    packages = [""],
    version = "0.0.3",
    description = "Sort, index and search BAM files by read ID",
    author = "David A. Parry",
    author_email = "david.parry@ed.ac.uk",
    url = "https://github.com/gantzgraf/index_bam_by_read_id",
    download_url = 'https://github.com/gantzgraf/index_bam_by_read_id/archive/0.0.3.tar.gz',
    license='MIT',
    test_suite='nose.collector',
    tests_require=['nose'],
    install_requires=[
          'pysam',
      ],
    classifiers = [
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
