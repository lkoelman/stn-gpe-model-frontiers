"""
A setuptools based setup module.

USAGE
-----

Build NMODL files:

    >>> python setup.py mechanisms

Install editable version (symlink to this directory):

    >>> python setup.py develop

Or, equivalently:

    >>> pip install -e path_to_repository
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# Custom build steps. You can either inherit from the base class
# distutils.cmd.Command or from one of the built-in classes:
# distutils.command.build<''/'_clib'/'_ext'/'_py'/'_scripts'>
from distutils.command.build import build as BuildCommand
# from distutils.cmd import Command as BuildCommand

# To use a consistent encoding
from codecs import open
from os import path
import os, platform
import subprocess

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    readme_contents = f.read()


BASEDIR = os.path.dirname(os.path.abspath(__file__))

ARCH = platform.machine()

class Build_NMODL(BuildCommand):
    """
    Try to compile NEURON NMODL mechanisms.
    """

    description = 'compile NEURON NMODL source files'

    # NOTE: can only have options if we subclass distutils.cmd.Command
    # user_options = [
    #     # The format is (long option, short option, description).
    #     ('nrnivmodl-path=', None, 'path to NEURON nrnivmodl executable'),
    # ]

    # def initialize_options(self):
    #     """
    #     Set default values for options.
    #     """
    #     BuildCommand.initialize_options(self)

    #     # Each user option must be listed here with their default value.
    #     self.nrnivmodl_path = ''

    #     # self.build_base = os.path.join(os.getcwd(), 'build-nmodl')
    #     # self.build_lib = os.path.join(os.getcwd(), 'build-nmodl')


    # def finalize_options(self):
    #     """
    #     Post-process options.
    #     """
    #     if self.nrnivmodl_path:
    #         assert os.path.isfile(self.nrnivmodl_path), (
    #             'Executable file {} does not exist.'.format(self.nrnivmodl_path))


    def run(self):
        """
        Run our custom build step
        """
        if getattr(self, 'nrnivmodl_path', False):
            nrnivmodl = self.nrnivmodl_path
        else:
            nrnivmodl = self._find_executable("nrnivmodl")

        if not nrnivmodl:
            print("Unable to find nrnivmodl. "
                  "You will have to compile NEURON .mod files manually.")
            return

        print("nrnivmodl found at", nrnivmodl)

        # Build mechanism files
        for root, dirs, files in os.walk(BASEDIR):
            if any((f.endswith('.mod') for f in files)) and not (
                root.endswith(ARCH)):
                # run `nrnivmodl` on our directory
                result, stdout = self._run_sys_command(nrnivmodl, root)

                if result != 0:
                    print("Unable to compile NMODL files in {dir}. Output was:\n"
                          "\t{output}".format(dir=root, output=stdout))
                else:
                    print("Successfully compiled NMODL files in {}.".format(root))



    def _find_executable(self, command):
        """
        Try to find an executable file.
        """
        path = os.environ.get("PATH", "").split(os.pathsep)
        cmd = ''
        for dir_name in path:
            abs_name = os.path.abspath(os.path.normpath(os.path.join(dir_name, command)))
            if os.path.isfile(abs_name):
                cmd = abs_name
                break
        return cmd


    def _run_sys_command(self, path, working_directory):
        p = subprocess.Popen(path, shell=True, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True,
                             close_fds=True, cwd=working_directory)
        result = p.wait()
        stdout = p.stdout.readlines()
        return result, stdout


# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    # This is the name of your project. The first time you publish this
    # package, this name will be registered for you. It will determine how
    # users can install this project, e.g.:
    #
    # $ pip install sampleproject
    #
    # And where it will live on PyPI: https://pypi.org/project/sampleproject/
    name='bgcellmodels',
    version='0.1.0',

    # This is a one-line description or tagline of what your project does. This
    # corresponds to the "Summary" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#summary
    description='A collection of Basal Ganglia cell and network models',  # Required

    # This is an optional longer description of your project that represents
    # the body of text which users will see when they visit PyPI.
    long_description=readme_contents,  # Optional
    # long_description_content_type='text/markdown',  # Optional (see note above)

    # This should be a valid link to your project's main homepage.
    url='https://bitbucket.org/lkmn_ucd/bg-cell-models',  # Optional

    # This should be your name or the name of the organization which owns the
    # project.
    author='Lucas Koelman',
    author_email='lucas.koelman@ucdconnect.ie',

    # Classifiers help users find your project by categorizing it.
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Neuroscience',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish
        'License :: OSI Approved :: GNU Lesser General Public '
        'License v3 (LGPLv3)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3.6',

        'Operating System :: POSIX',
        'Natural Language :: English',
    ],

    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    keywords='computational neuroscience simulation neuron',  # Optional

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    packages=['bgcellmodels'],
    # packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required

    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    install_requires=[
        'numpy',
        'matplotlib',
        'pint',
        'transforms3d',
        'plyfile',
        'lazyarray',
        'pint',
        # Manual install specific versions:
        # 'PyNN',
        # 'BluePyOpt',
        # 'elephant',
    ],

    # List additional groups of dependencies here (e.g. development
    # dependencies). Users will be able to install these using the "extras"
    # syntax, for example:
    #
    #   $ pip install sampleproject[dev]
    #
    extras_require={  # Optional
        'optimisation': ['cython', 'numba', 'PySpike'],
        'devtools': ['nb_conda', 'nbstripout', 'jupyter_contrib_nbextensions'],
    },

    # If there are data files included in your packages that need to be
    # installed, specify them here.
    #
    # If using Python 2.6 or earlier, then these have to be included in
    # MANIFEST.in as well.
    package_data={  # Optional
        'bgcellmodels': [
            'bgcellmodels/models/STN/GilliesWillshaw/sth-data/*',
            'bgcellmodels/models/GPe/Gunay2008/morphology/*'
        ],
    },

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    #
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],  # Optional


    # Extra build steps defined as subclasses of distutils.cmd.Command:
    # These are invoked as `python setup.py <key>`
    cmdclass={'mechanisms': Build_NMODL},

    # List additional URLs that are relevant to your project as a dict.
    #
    # This field corresponds to the "Project-URL" metadata fields:
    # https://packaging.python.org/specifications/core-metadata/#project-url-multiple-use
    #
    # Examples listed include a pattern for specifying where the package tracks
    # issues, where the source is hosted, where to say thanks to the package
    # maintainers, and where to support the project financially. The key is
    # what's used to render the link text on PyPI.
    project_urls={  # Optional
        'Lab Website': 'https://www.neuromuscularsystemsucd.info/',
        'University Website': 'http://www.ucd.ie/biomedicalengineering/research/neuroengineeringoptics/neuromuscularsystemsandneuralengineeringlab/',
        'Personal Website': 'https://lkoelman.github.io/',
        'Source': 'https://bitbucket.org/lkmn_ucd/bg-cell-models/src',
        # 'Bug Reports': 'https://github.com/pypa/sampleproject/issues',
    },
)
