./bin/yscenegen tests/01_terrain/terrain.json -o outs/01_terrain/terrain.json --terrain object
./bin/yscenegen tests/02_displacement/displacement.json -o outs/02_displacement/displacement.json --displacement object
./bin/yscenegen tests/03_hair1/hair1.json -o outs/03_hair1/hair1.json --hairbase object --hair hair
./bin/yscenegen tests/03_hair2/hair2.json -o outs/03_hair2/hair2.json --hairbase object --hair hair --hairlen 0.005 --hairstr 0
./bin/yscenegen tests/03_hair3/hair3.json -o outs/03_hair3/hair3.json --hairbase object --hair hair --hairlen 0.005 --hairstr 0.01
./bin/yscenegen tests/03_hair4/hair4.json -o outs/03_hair4/hair4.json --hairbase object --hair hair --hairlen 0.02 --hairstr 0.001 --hairgrav 0.0005 --hairstep 8
./bin/yscenegen tests/04_grass/grass.json -o outs/04_grass/grass.json --grassbase object --grass grass

./bin/yscenetrace outs/01_terrain/terrain.json -o out/01_terrain.jpg -s 256 -r 720
./bin/yscenetrace outs/02_displacement/displacement.json -o out/02_displacement.jpg -s 256 -r 720
./bin/yscenetrace outs/03_hair1/hair1.json -o out/03_hair1.jpg -s 256 -r 720
./bin/yscenetrace outs/03_hair2/hair2.json -o out/03_hair2.jpg -s 256 -r 720
./bin/yscenetrace outs/03_hair3/hair3.json -o out/03_hair3.jpg -s 256 -r 720
./bin/yscenetrace outs/03_hair4/hair4.json -o out/03_hair4.jpg -s 256 -r 720
./bin/yscenetrace outs/04_grass/grass.json -o out/04_grass.jpg -s 256 -r 720 -b 128

######################
./bin/yscenegen tests/extra1/extra1.json -o outs/extra1/extra1.json --extra1 object

./bin/yscenetrace outs/extra1/extra1.json -o out/uniformed_points.jpg -s 256 -r 720

./bin/yscenegen tests/extra2/extra2.json -o outs/extra2/extra2.json --extra2 object

./bin/yscenetrace outs/extra2/extra2.json -o out/lunar_surface.jpg -s 256 -r 720
