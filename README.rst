===========================
Proposal for Video Compression
===========================

Basic information
=================

`My GitHub repository (URL) <https://github.com/jspss93094/nsd_TermProject>`

Use C++ with Python to implement the video compression, and the file type can be various.

Problem to solve
================

Frequently, the videos are taken up much space if we wanted to save it, download it, or do something.
It must took a long time to deal with a video. Even its size will be larger when the resolution is higher.

Some problems I met when I want to compress the video:

1. It often needed to install the software.
2. It took a long time to compress, but the result was not better than before.
3. The resolution became too low to watch.

Perspective users
=================

The person who :

1. Wanted to compress the large video.
2. Do not want to install the software.
3. Can not upload the file online because of the large size.
4. Found the results of compression is not well. 

System architecture
===================

It also have many modules can use, like FFMPEG. But I want to realize the process and to implement it.
If some process in Python is better than in C++, I will embed some Python.

I just use C++ to create the function to complie and execute, or the file which can be imported in Python.

The function to implement : 

1. Get the input.
2. Decode the input video.
3. Do the calculation of compression.
4. Encode it to be the video type.
5. Output the result and save.

API description
===============

If users want to use my program in the program by themself, it can be easy to 

include the file in C++. Just download my file and include it in your header.

I will also provide the file which can be imported in Python. It can use by 

just call the functions.

Maybe I will make the graphical user interface which is easier to use by the general users.

Engineering infrastructure
==========================

Use make to compile the cpp, it will genarate the files :

1. Import module in Python.
2. The executable file only show in the console window.
3. The executable file for general users, it will show the graphical user interface.

Version control :

1. When the function is updated to the certain degree, I will upload the new version both 
   to the google drive and github branch.
2. When it updated, the contents of updating details will show in the update.txt file.
3. In order to control the version branch, I will let all versions of the same function save 
   in the different directory.
4. The github branch will just show the function name.
5. The main branch will always be updated for the version which can be executed definitely.

Schedule
========

- week 7 (11/01): video file code type
- week 8 (11/08): midterm
- week 9 (11/15): get input function
- week 10(11/22): calculation implement
- week 11(11/29): calculation implement
- week 12(12/06): calculation implement
- week 13(12/13): output and create module
- week 14(12/20): graphical user interface
- week 15(12/27): prepare presentation

References
==========

`An Introduction to Video Compression in C/C++ 
<https://www.yumpu.com/en/document/read/19876586/an-introduction-to-video-compression-in-c-c->`

`Learning for Video Compression <https://arxiv.org/pdf/1804.09869.pdf>`

`Library <https://www.itread01.com/content/1547186605.html>`
