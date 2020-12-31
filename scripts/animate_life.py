import glob
from PIL import Image
from functools import partial

# read in metadata file entries to a list of lists
print("Loading simulation metadata...")
fname = "./outfiles/metadata.txt"
file = open(fname, "r")
metadata = []
for line in file:
    line = line.split()
    metadata.append(line)

# parse metadata variables
rows = int(metadata[0][1])
cols = int(metadata[0][2])
num_procs = int(metadata[1][1])
iprocs = int(metadata[2][1])
jprocs = int(metadata[2][2])
num_gen = int(metadata[3][1])
is_periodic = int(metadata[4][1])
freq_dump = int(metadata[5][1])


def create_animation():
    palette = [0, 0, 0, 255, 255, 255]
    palette = palette + [0]*(768-len(palette))

    paths = glob.glob('./outfiles/dump*')
    imgs = []
    for path in sorted(paths, key=lambda x: int(x.split("p")[1])):
        with open(path, 'rb') as ifile:
            for data in iter(partial(ifile.read, rows*cols), b''):
                img = Image.frombuffer('L', (cols, rows), data)
                img.putpalette(palette)
                imgs.append(img)
        imgs[0].save('animation.gif', save_all=True,
                     append_images=imgs[1:], loop=0)


# write and save the animation
print("Writing ./outfiles/animation.gif...")
create_animation()
print("\nAnimation ready.")
