# Ray Tracer

Rust implementation for [The Ray Tracer Challenge](http://raytracerchallenge.com/).

## Usage

```text
The Ray Tracer Challenge CLI

Usage: raytracer [OPTIONS] --scene <FILE>

Options:
      --scene <FILE>  A yaml description of the scene to render
      --obj <FILE>    Optional obj models to add to the scene
      --out <OUT>     Optional output ppm file, defaults to stdout
      --ppm <FILE>    Optional ppm textures to use as material
  -h, --help
```

Note that multiple `obj` models can be added/combined in a scene.

## Examples

Render a basic scene, output it to stdout, pipe it to [ImageMagick](https://imagemagick.org/):

```
raytracer --scene samples/scenes/basic_scene.yaml | magick display
```

Render the book's cover and output it to a `.ppm` file:

```
raytracer --scene samples/scenes/cover.yaml --out cover.ppm
```

Render an `obj` model and output it to a `.ppm` file:

```
raytracer --scene samples/scenes/space_ship.yaml --obj samples/obj/space_ship.obj --out space_ship.ppm
```

Render a teapot over a space ship !

```
raytracer --scene samples/scenes/space_teapot.yaml \
          --obj samples/obj/space_ship.obj \
          --obj samples/obj/teapot_low.obj \
          --out space_teapot.ppm
```
