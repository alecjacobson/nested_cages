Small program for view results of the "Nested Cages" SIGGRAPH Submission.

Mac OS X 10.10 only

Double-click on the NestView.app to launch the viewer on a default set of cages
around the bunny. You'll see the bunny loaded and then the program will pause
as it creates the "clip surface" between layers. Once finished you can:

[Click] and [drag]             Rotate scene.
[Click] and [drag] on widget   Rotate clip plane.
[Scroll]                       Zoom in and out.
[⇧ Scroll]                     Push clip plane in and out.
[Z]                            Snap clip plane rotation to canonical view.
[z]                            Snap rotation to canonical view.
[⌘ Z]                          Undo.
[⇧ ⌘ Z]                        Redo.
[^C,ESC]                       Exit.
[-]/[+]                        Zoom out/zoom in.
[h]                            Toggle whether to clip finest model.
[O]                            Toggle visibility of rotation widget.
[o]                            Load model from single .obj out of sequence
[.]/[,]                        Decrease/increase 'staggering' of clipping plane.
[<]/[>]                        Decrease/increase clipping plane offset.
[(]/[)]                        Decrease/increase number of layers visible.
[l]                            Toggle visibility of wireframe lines.
[S,s]                          Save camera and cut meshes.
[SPACE]                        Recompute clipping planes.

Load a new model by pressing 'o' or clicking on 'load .obj'. Then select one of
the .obj's from the sequences with filenames of the form `path/to/model_0.obj`,
`path/to/model_1.obj`,...,`path/to/model_[n].obj`. Only .objs will be parse
(though all files are allowed by the dialog).

The program may appear to hang while computing the clip surface. Please be
patient.
