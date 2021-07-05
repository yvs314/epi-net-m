
## Contributing
1. The `master` branch **must** be kept in functioning state at all times. 

2. `git branch` yourself a private work space, and, when the code is done and tested, make a **pull request** and tag Yaroslav Salii @yvs314 to review and merge it into `master`.

## Code and *commit* style
1. Prefer meaningful commit messages. Say, “Updated script.m” is worse than “Fix typos in comments”. Preferably start the messages with a verb in neutral form, e.g. add, fix, upd, rehaul, etc. 

Add an *extended description* whenever you fix a notable bug. Write up an **issue** if there is a problem.

2. Please **do not** leave the old versions of functions in the code. If you feel like retaining commented-out code, consider _stashing_ it somewhere or keeping it in a **bit bucket** after the end of working code.

3. Isolate *source code* from *input* and *output*.

4. Use _Git Large File Storage_ for data, figures, and anything else that is not code. See the [this page](https://git-lfs.github.com/) for how to get started on it.

Randall Munroe's [take](https://xkcd.com/1296/) on commit messages:

![uts](https://imgs.xkcd.com/comics/git_commit.png)
