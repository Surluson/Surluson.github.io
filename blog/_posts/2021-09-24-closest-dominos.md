---
layout: post
title: A home delivery optimization problem
truncated_preview: true
excerpt_separator: <!--more-->
usemathjax: true
published: true
---

![apple_pie](/blog/visualize_svd_extra/apple_pie.png "A "great" looking apple pie")

<div class="message">
    Recently I published a paper utilizing the Singular Value
    Decomposition (SVD), a well-known matrix decomposition method.
    With my first blog post on my site, I wish to explain how to
    intuitively visualize the results the singular value decomposition
    provides using a simple example.
</div>

<!--more-->

## Motivation

During one of my lazy evenings I didn't feel like cooking so I ordered a Domino's pizza which would be delivered to my doorsteps.
Being a moderately impatient person, I started wondering how far the closest Domino's location was from my home relative to other places in Reykjavik.
The home delivery involved driving into my district, which requires the driver to drive further away to enter the district and essentially turn back to get to my house, so the naive method of comparing the Euclidean distance (or the aerial distance if you will) between the Domino's location and my home would probably be insufficient.
A better way would be to take the road network in Reykjavik into consideration, and figure out the shortest path from point A (my home) and point B (Domino's) and measure the distance on that path.

The problem of finding the shortest path to a Domino's for myself would be a pretty boring task, as I have a pretty good notion of the shortest path to the closest Domino's, having been driving in my district for years.
Gathering the shortest distance to any Domino's for every single home in Reykjavik would be a more interesting task, because it would require me to (1) measure the distance between points A' (all homes in Reykjavik) and points B' (all Domino's locations in Reykjavik), and (2) figure out the path with the shortest distance to find the closest Domino's and compare that information with other homes in Reykjavik.

Fortunately, there exists an open-source project which contains essentially all the data I need for this project: [OpenStreetMap](https://www.openstreetmap.org/about)

## Geospatial Data and Working with Open Street Map

Open Street Map (OSM) is a gigantic open-source project founded in 2004 and is a global database of geographic locations, roads, landmarks and more.
Being open-source, users from all over are work on updating the data, and in some cases it means that the data in some places might be underrepresented in the database.

In my case, I got a `.osm` file of Iceland (a whopping gigabyte of data), but to keep my computer from being fried I reduced the scope of the geodata to only view Reykjavik, which reduced the size of the file down to ~100 MB, with 800431 distinct points, and 103230 roads mapped (actually it turned out to be 103219 roads after removing unconnected roads that popped up after reducing the scope of the data).

Another issue I ran into was that the location and data for the Domino's restaurants was not up-to-date in some cases.
One restaurant was completely missing from the data while one restaurant which has been closed down for years was still present in the data (which I've since updated so you won't see these issues).

To work with the geodata I used a library called [OpenStreetMapX.jl](https://github.com/pszufe/OpenStreetMapX.jl) in Julia, which converts the geodata into a graph that I can exploit using various algorithms.

## The Reykjavik Graph

The capital area consists of 7 municipals, but for this project I will ignore one completely (Kjósarhreppur) and another partially (Mosfellsbær).
The reason being that Kjósarhreppur exists far away from what people would regularly call the city limits, and there is not an interesting road network compared to the densely populated area that I will consider.
The part of Mosfellsbær that I omitted was mostly an error on my part, but when I shrunk the scope of the geodata I accidentally skipped one neighborhood, but simple reasoning can be used to extrapolate the information I gather to that specific neighborhood.

After updating the geodata I needed, the map contained 18 Domino's locations.
One issue I ran into was that the Domino's locations were not usually a part of the road network, but rather "shopping" nodes the roads lead to.
Without increasing the scope of my project significantly, I got rid of all nodes and edges that did not belong to the road network, but prior to that I located the road network nodes that were closest to the Domino's locations and marked them.
That way I was able to get rid of a significant portion of the nodes while still keeping a good estimate of where the Domino's locations were.

<div>
    <iframe src="https://surluson.github.io/map.html" height="560" width="560" allowfullscreen="" frameborder="0">
    </iframe>
</div>