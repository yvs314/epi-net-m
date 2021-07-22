
## Contributing
1. The `master` branch **must** be kept in functioning state at all times. 

2. Whenever adding new features, `git branch` yourself a private work space, and, when the code is done and tested, make a **pull request** and tag the maintainer, Yaroslav Salii @yvs314, to review and merge it into `master`.

3. It is OK to commit directly to `master` when updating the **documentation**, which includes both separate .md files and _comments_.

## Code and *commit* style
1. Prefer meaningful commit messages. Say, “Updated script.m” is worse than “Fix typos in comments”. Preferably start the messages with a verb in neutral form, e.g. add, fix, upd, rehaul, etc. 

Add an *extended description* whenever you fix a notable bug. Write up an **issue** if there is a problem.

2. Prefer **not** to retain the old/experimental code inline, commented-out. If you feel like retaining commented-out code, consider _stashing_ it somewhere or keeping it in a **bit bucket** after the end of working code.

3. Isolate *source code* from *input* and *output*.

4. Use _Git Large File Storage_ for data, figures, and anything else that is not code. See the [this page](https://git-lfs.github.com/) for how to get started on it.

Randall Munroe's [take](https://xkcd.com/1296/) on commit messages:

![uts](https://imgs.xkcd.com/comics/git_commit.png)

### Semantic Versioning

Use version format  v.<major>.<minor>\[.<fix>\], where **major** starts at 0, increases to 1 when feature-complete and more on major rewrites, **minor** tracks complete features, and **fix** is optional, only appearing on squashing a bug or during eureka moments that don't qualify as a complete feature. 

It is fine to use hexadecimals `a`–`f` to stretch single-figure numbers from ten `0`–`9` to sixteen `0`–`f`.

### Inline Docs and Attribution
Start each code file with authorship `Author: <name>, <year>`. Then describe the purpose of this file in a couple of sentences, and, if applicable, add invocation instructions or sample as well.

### Changelog

Reflect _noteworthy_ changes at the start of the file, directly below the _Inline Docs and Attribution_. Use _Semantic Versioning_ as described above. Suggested format is ISO-style YYYY-MM-DD date, version number, version description, separated by tabs or spaces, e.g.
```
Oboe.jl
2020-12-18   v.0.7: wrote pairs-to-matrix xform by hand (vs. `unstack`, which was unpredictable)
2020-12-30   v.0.8: added a beta node-to-node daily air passenger computation
2021-01-03   v.0.8.1: a half-baked Main(), look in bit bucket. Tested on 3K by-county!
```
Version descriptions can be slightly longer than commit messages.
