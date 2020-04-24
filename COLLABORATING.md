Lemme describe (1) code and commit style and (2) how to properly push your work here.

## Code and *commit* style
1. Please write meaningful commit messages, say, “Updated script.m” is worse than “Fixed a few typos in comments”. 
If you fix a bug, an *extended description* is most welcome, to make your fix more specific.
Here's what Randall Munroe's [take](https://xkcd.com/1296/) on it:
![uts](https://imgs.xkcd.com/comics/git_commit.png)


2. Please **do not** leave the code you modified as a comment. With `git`'s built-in changes tracking it is *redundant* and actually makes it *harder* to compare. If you want to compare different versions, see the history at `GitHub`'s web interface (or jockey with `git diff` in your console).

## Adding your work
1. The `master` branch **must** be kept in functioning state at all times. 

That is, whenever you push to `master`, you must test that no functionality is broken, even if you did not touch it directly. Say, I modify the *initial values' format* and update the input subroutines. I must still verify that the *numeric solution* runs fine—just in case I have inadvertently renamed the wrong variable or something.

2. For *major* changes, do not commit directly to `master`.
Make yourself a new `branch`, iron out the kinks in it at your leisure, and then file a `pull request` to merge it into the `master`.

I also think it's useful to review each other's code, that is, make a `pull request` but do not `merge` it yourself, ask your colleague to do it. But we'll see if it will work.
