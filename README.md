# D-geodiversity

This blogram finds the D-geodiversity of the two pathways.
For a detailed explanation of D-geodiversity, please see Ref.
The usage of this program is Enter the node of interest in node_init.txt. 
Then, by executing main.cpp, the D-geodiversity of the node and any two points of interest will all be stored in distance.txt

The topology used in this project was created for 50 German cities.
Node_info.txt contains latitude and longitude information of 50 German cities.
In Links_info.txt, the first column contains the link numbers, the second and third columns contain the connected node numbers, and the fourth column contains the length of the links.

Note: The input coordinates must be latitude and longitude because this program calculates the distance from latitude and longitude.

Reference:
D. Santos, T. Gomes, and D. Tipper, “SDN Controller Placement With Availability Upgrade Under Delay and Geodiversity Constraints,” IEEE Transactions on Network and Service Management, vol. 18, no. 1, pp. 301–304, Mar. 2021
