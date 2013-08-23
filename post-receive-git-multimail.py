#! /usr/bin/python2.7

"""Example post-receive hook based on git-multimail.

This script is a simple example of a post-receive hook implemented
using git_multimail.py as a Python module.  It is intended to be
customized before use; see the comments in the script to help you get
started.

It is possible to use git_multimail.py itself as a post-receive or
update hook, configured via git config settings and/or command-line
parameters.  But for more flexibility, it can also be imported as a
Python module by a custom post-receive script as done here.  The
latter has the following advantages:

* The tool's behavior can be customized using arbitrary Python code,
  without having to edit git_multimail.py.

* Configuration settings can be read from other sources; for example,
  user names and email addresses could be read from LDAP or from a
  database.  Or the settings can even be hardcoded in the importing
  Python script, if this is preferred.

This script is a very basic example of how to use git_multimail.py as
a module.  The comments below explain some of the points at which the
script's behavior could be changed or customized.

"""

import sys
import os
import re
import optparse

# If necessary, add the path to the directory containing
# git_multimail.py to the Python path as follows.  (This is not
# necessary if git_multimail.py is in the same directory as this
# script):

LIBDIR = '/var/git/cgal.git/hooks/git-multimail/git-multimail'
#LIBDIR = '/home/lrineau/Git/CGAL-hooks-on-server/git-multimail/git-multimail'
sys.path.insert(0, LIBDIR)

import git_multimail


# It is possible to modify the output templates here; e.g.:

#git_multimail.FOOTER_TEMPLATE = """\
#
#-- \n\
#This email was generated by the wonderful git-multimail tool.
#"""

git_multimail.REF_CREATED_SUBJECT_TEMPLATE = (
    '%(emailprefix)s%(refname_type)s %(short_refname)s created.'
)
git_multimail.REF_UPDATED_SUBJECT_TEMPLATE = (
    '%(emailprefix)s%(refname_type)s \'%(short_refname)s\' updated. %(packages_short)s'
)
git_multimail.REF_DELETED_SUBJECT_TEMPLATE = (
    '%(emailprefix)s%(refname_type)s %(short_refname)s deleted.'
)

git_multimail.REFCHANGE_INTRO_TEMPLATE = """\
repo:     %(repo_shortname)s
reftype:  %(refname_type)s
refname:  %(short_refname)s
pusher:   %(pusher)s
packages: %(packages)s

https://scm.cgal.org/gitweb/?p=%(repo_shortname)s.git;a=shortlog;h=%(newrev)s

"""

git_multimail.REVISION_HEADER_TEMPLATE = """\
To: %(recipients)s
Subject: %(emailprefix)s%(num)02d/%(tot)02d: %(oneline)s %(packages_short)s
MIME-Version: 1.0
Content-Type: text/plain; charset=%(charset)s
Content-Transfer-Encoding: 8bit
From: %(fromaddr)s
Reply-To: %(reply_to)s
In-Reply-To: %(reply_to_msgid)s
References: %(reply_to_msgid)s
X-Git-Repo: %(repo_shortname)s
X-Git-Refname: %(refname)s
X-Git-Reftype: %(refname_type)s
X-Git-Rev: %(rev)s
X-Git-CGALpackages: %(packages)s
Auto-Submitted: auto-generated
"""

git_multimail.REVISION_INTRO_TEMPLATE = """\
repo:     %(repo_shortname)s
reftype:  %(refname_type)s
refname:  %(short_refname)s
pusher:   %(pusher)s
packages: %(packages)s

https://scm.cgal.org/gitweb/?p=%(repo_shortname)s.git;a=commitdiff;h=%(rev)s

"""

# The template used in summary tables.  It looks best if this uses the
# same alignment as TAG_CREATED_TEMPLATE and TAG_UPDATED_TEMPLATE.
# TODO-FUTURE packages=packages_outb(shorten_packages(get_packages(self.new))),
#git_multimail.BRIEF_SUMMARY_TEMPLATE = """\
#%(action)10s  %(rev_short)-9s %(text)s %(packages_short)s
#"""

git_multimail.FOOTER_TEMPLATE = ( '' )
git_multimail.REVISION_FOOTER_TEMPLATE = ( '' )


# Specify which "git config" section contains the configuration for
# git-multimail:
#config = git_multimail.Config('multimailhook')

# Select the type of environment:
#environment = git_multimail.GenericEnvironment(config)
#environment = git_multimail.GitoliteEnvironment(config)

def get_packages(rev):
    packages = list(git_multimail.read_git_lines(['diff-tree', '--no-commit-id', '--name-only', '-r', str(rev)], keepends=False))
    # each file in root dir (i.e. without "/") will be replaced by "<root>"
    packages = [re.sub(r'^(?!.*\/).*$', r'<root>', p) for p in packages]
    # next strip off path after first "/"
    packages = [re.sub(r'\/.*', r'', p) for p in packages]
    return sorted(set(packages))

def get_packages2(rev1, rev2):
    # question: do we need  ['%s..%s' % (rev1, rev2,)] instead?
    packages = list(git_multimail.read_git_lines(['diff-tree', '--name-only', '-r', rev1, rev2], keepends=False))
    # each file in root dir (i.e. without "/") will be replaced by "<root>"
    packages = [re.sub(r'^(?!.*\/).*$', r'<root>', p) for p in packages]
    # next strip off path after first "/"
    packages = [re.sub(r'\/.*', r'', p) for p in packages]
    return sorted(set(packages))

def shorten_packages(packages):
    return packages if len(packages) < 4 else ["--multiple--"]

def packages_out(packages):
    return ", ".join(packages)

def packages_outb(packages):
    return "{" + ", ".join(packages) + "}"


class CgalScmEnvironment(
        git_multimail.ProjectdescEnvironmentMixin,
        git_multimail.ConfigMaxlinesEnvironmentMixin,
        git_multimail.ConfigFilterLinesEnvironmentMixin,
        git_multimail.PusherDomainEnvironmentMixin,
        git_multimail.ConfigOptionsEnvironmentMixin,
        git_multimail.GenericEnvironmentMixin,
        git_multimail.ConfigRecipientsEnvironmentMixin,
        git_multimail.GitoliteEnvironmentMixin,
        git_multimail.Environment,
):
    pass

    def get_pusher(self):
        return os.environ.get('GIT_NAME', 'Unknown user')

    def get_pusher_email(self):
        return self.get_pusher() + (' <%s@users.gforge.inria.fr>' % os.environ.get('GIT_USER', 'unknown_user'))

    def get_fromaddr(self):
        return self.get_pusher_email()

    def get_sender(self):
        return 'git@scm.cgal.org'

    def get_refchange_recipients(self, refchange):
        """Abused to set packages for refchange summary email"""

        # need to init here already, as called very early in git-multimail
        self.init_values()

        old = refchange.old
        new = refchange.new

        if old.sha1 == None or new.sha1 == None:
            self._values['packages'] = "n/a"
            self._values['packages_short'] = ""
        else:
            packages = get_packages2(str(old), str(new))
            self._values['packages'] = packages_out(packages)
            self._values['packages_short'] = packages_outb(shorten_packages(packages))

        return self._get_recipients(self.config, 'refchangelist', 'mailinglist')

    def get_revision_recipients(self, revision):
        """Abused to set packages for commit email of revision"""

        packages = get_packages(revision.rev)
        self.init_values()
        self._values['packages'] = packages_out(packages)
        self._values['packages_short'] = packages_outb(shorten_packages(packages))
        return self._get_recipients(self.config, 'refchangelist', 'mailinglist')

    def init_values(self):
        if self._values is None:
            values = {}

            for key in self.COMPUTED_KEYS:
                value = getattr(self, 'get_%s' % (key,))()
                if value is not None:
                    values[key] = value

                #values['packages'] = ''
                #values['packages_short'] = ''

            self._values = values

    def get_values(self):
        """Return a dictionary {keyword : expansion} for this Environment.

        This method is called by Change._compute_values().  The keys
        in the returned dictionary are available to be used in any of
        the templates.  The dictionary is created by calling
        self.get_NAME() for each of the attributes named in
        COMPUTED_KEYS and recording those that do not return None.
        The return value is always a new dictionary."""

        self.init_values()

        return self._values.copy()


def main(args):
    parser = optparse.OptionParser(
        description=__doc__,
        usage='%prog [OPTIONS]\n   or: %prog [OPTIONS] REFNAME OLDREV NEWREV',
        )

    parser.add_option(
        '--stdout', action='store_true', default=False,
        help='Output emails to stdout rather than sending them.',
        )
    parser.add_option(
        '--recipients', action='store', default=None,
        help='Set list of email recipients for all types of emails.',
        )

    (options, args) = parser.parse_args(args)

    config = git_multimail.Config('multimailhook')
    environment = CgalScmEnvironment(config=config)

    if not config.get_bool('enabled', default=True):
        sys.exit(0)

    try:
        if options.stdout:
            mailer = git_multimail.OutputMailer(sys.stdout)
        else:
            mailer = git_multimail.choose_mailer(config, environment)


        # Dual mode: if arguments were specified on the command line, run
        # like an update hook; otherwise, run as a post-receive hook.
        if args:
            if len(args) != 3:
                parser.error('Need zero or three non-option arguments')
            (refname, oldrev, newrev) = args
            git_multimail.run_as_update_hook(environment, mailer, refname, oldrev, newrev)
        else:
            git_multimail.run_as_post_receive_hook(environment, mailer)
    except git_multimail.ConfigurationException, e:
        sys.exit(str(e))

if __name__ == '__main__':
    main(sys.argv[1:])
