#include <iostream>
#include <algorithm>
#include <vector>
#include <ctime>
#include <stdbool.h>
#include <cmath>
#include <time.h>
#include <fstream>

using namespace std;

#define LINKER 88
#define NODER 50
#define D 50
#define MIN_DISTANCE 10000


typedef struct{
	double lng;
    double lat;
}POINT;


typedef struct{
	int id;
	int a;
	int b;
	int distance;
}LINK;



int Dijkstra(int current, vector<vector<int>> &short_node, int *route, int network[][NODER][2]){
	int pDist[NODER];/* Array to set the shortest distance from the starting point to each location */
	int pRoute[NODER];
	bool pFixed[NODER];/* Array to identify if the shortest distance from the starting point to each point has been determined */
	int sPoint,i,j,newDist,shortPoint = 0;
	int max_distance = 0;

	vector<vector<int>> data(NODER, vector<int> (2));

	/* Stores the initial value for the shortest distance from the starting point to the destination (no need to change) */

	for(i=0;i<NODER;i++){
		//sRoute[current][i] = -1; /* Stores the initial value in the point number of the point on the shortest route. */
		pDist[i]=99999; /* Stores an initial value for the shortest distance from the starting point to each point. */
		pFixed[i]=false; /* Stores the initial value in the fixed state of the shortest distance at each location. */
	}

	pDist[current]=0;/* Set 0 to the shortest distance from the starting point to the starting point itself */

	while(true){ /* shortest path search process */
		i=0;
		while(i<NODER){/* Find one undetermined point. */
			if(pFixed[i]==0){
				break; /* Exit from re-interior repetition */
			}
			i=i+1;
		}

		if(i==NODER){ /* If the shortest route from the starting point to all points is determined. */
			break; /* Exit the shortest route search process. */
		}

		for(j=i+1;j<NODER;j++){ /* Find the point where the shortest distance is shorter */
			if((pFixed[j]==0) && (pDist[j] < pDist[i])){
				i=j;
			}
		}

		sPoint=i;
		pFixed[sPoint]=true; /* Determine the shortest distance from the point of departure */
		data.at(shortPoint).at(0) = sPoint;
		data.at(shortPoint).at(1) = pDist[sPoint];
		shortPoint++;
		
		
		for(j=0;j<NODER;j++){
			if((network[sPoint][j][1]>0) && (pFixed[j]==0)){
				newDist=pDist[sPoint]+network[sPoint][j][1];
				if(newDist<pDist[j]){
					pDist[j]=newDist;
					route[j]=network[sPoint][j][0];
				}
			}
		}
	}
	sort(data.begin(),data.end(),[](const vector<int> &alpha,const vector<int> &beta){return alpha[1] < beta[1];});
	short_node = data;
	max_distance = short_node.back().back();
	return max_distance;
}

int Distance(POINT a, POINT b) {

	double lat1 = a.lat;
	double lng1 = a.lng;
	double lat2 = b.lat;
	double lng2 = b.lng;

    // Pi (3.1415926..)
    const double pi = 3.14159265359;

    // Convert latitude and longitude to radians
    double rlat1 = lat1 * pi / 180;
    double rlng1 = lng1 * pi / 180;
    double rlat2 = lat2 * pi / 180;
    double rlng2 = lng2 * pi / 180;


    // Find the central angles (radians) of two points
    double c =
      sin(rlat1) * sin(rlat2) +
      cos(rlat1) * cos(rlat2) *
      cos(rlng1 - rlng2);
    double rr = acos(c);

    // Earth's equatorial radius (meter)
    const double earth_radius = 6378140;

    //  Distance between two points (meter)
    double dis = round((earth_radius * rr) / 1000) ;


	int distance = dis;
    return distance;
}

POINT intersection(POINT a, POINT b, POINT c){
	POINT inter;
	double x,y;//point of intersection
	double d = (a.lat - b.lat) / (a.lng - b.lng);//slope
	double e = (-a.lng) * d + a.lat;//section


	x = (d * (c.lat - e) + c.lng) / (pow(d,2) + 1);
	y = d * ((d * (c.lat - e) + c.lng) / (pow(d,2) + 1)) + e;

	inter.lng = x;
	inter.lat = y;

	double dis_ab = Distance(a,b);
	double dis_ac = Distance(a,inter);
	double dis_bc = Distance(b,inter);

	if(dis_ab < dis_ac){
		return b;
	}else if(dis_ab < dis_bc){
		return a;
	}else{
		return inter;
	}

}

int min_geodiversity(int i, int route[][NODER],LINK *link, POINT *point,int node_point1, int node_point2){
	int a,b,link_point1,link_point2,distance,min_distance = MIN_DISTANCE;
	int flag;
	
	POINT inter;

	while(node_point1 != i){
		link_point1 = route[i][node_point1];//Stores link numbers
		while(node_point2 != i){
			link_point2 = route[i][node_point2];//Stores link numbers
			a = link[link_point2].a;
			b = link[link_point2].b;//Stores the node number connected to link_point2
			// Set the distance to zero if the same route is taken
			if(link_point1 == link_point2 || a == node_point1 || b == node_point1){
				return 0;
			}

			
			inter = intersection(point[a],point[b],point[node_point1]); //Calculation of intersections
			distance = Distance(inter,point[node_point1]); //Distance Calculation
			if(distance < min_distance){
				min_distance = distance;
			}

			if(link[link_point2].a == node_point2){
				node_point2 = link[link_point2].b;
			}else{
				node_point2 = link[link_point2].a;
			}
		
		}


		if(link[link_point1].a == node_point1){
			node_point1 = link[link_point1].b;
		}else{
			node_point1 = link[link_point1].a;
		}


	}
	return min_distance;
}

int main(void){
	char fname1[] = "Links_info.txt";
	char fname2[] = "Node_info.txt";
	char fname3[] = "init_node.txt";
	char ch[50];
	
	int i,j,k;
	int distance;
	int network[NODER][NODER][2],route[NODER][NODER];
	double n,m;

	int i1,i2,p,sedai = 0,x = 0,nc = 0,penalty = 0,penalty_cc = 0,penalty_link = 0, con_ok,link_point,link_point1,link_point2,node_point,node_point1,node_point2,a,b,distance_link,distance_node,distance_nodea,distance_nodeb,flag = 0,max_distance = 0,max_dis[NODER],min_distance = 10000,min_distance1,min_distance2;
	int cross_point,delay,min_point,counter = 0;
	int node[NODER];

	FILE *fp1,*fp2,*fp3;
	vector<vector<vector<int>>> short_node(NODER,vector<vector<int>>(NODER, vector<int> (2)));
	vector<int> primary_c(2);
	
	LINK link[LINKER];//Link Information
  	POINT point[NODER],inter;

	ofstream outputfile("distance.txt");

	srand((unsigned)time(NULL));
    
	fp1 = fopen(fname1, "r");
	fp2 = fopen(fname2, "r");
	fp3 = fopen(fname3, "r");
	

	if(fp1 == NULL) {
		cout << fname1 << "file not open!" << endl;
		return -1;
	}

	if(fp2 == NULL) {
		cout << fname2 << "file not open!" << endl;
		return -1;
	}

	if(fp3 == NULL) {
		cout << fname3 << "file not open!" << endl;
		return -1;
	}

   

  	i = 0;
	for(i = 0; i < NODER; i++){
		node[i] = 0;
		for(j = i; j < NODER ; j++){
			for(int k = 0; k < 2; k++){
				if(i == j){
					network[i][i][k] = 0;
				}else{
					network[i][j][k] = -1;
					network[j][i][k] = -1;
				}
			}
		}
	}
	i = 0;
	while(fscanf(fp1,"%d %d %d %d",&link[i].id,&link[i].a,&link[i].b,&link[i].distance) != EOF){
		link[i].id--;
		link[i].a--;
		link[i].b--;
		network[link[i].a][link[i].b][0] = link[i].id;
		network[link[i].b][link[i].a][0] = link[i].id;
		network[link[i].a][link[i].b][1] = link[i].distance;
		network[link[i].b][link[i].a][1] = link[i].distance;
		i++;
	}

    i = 0;

    while (fscanf(fp2,"%d %s %lf %lf",&i,ch,&n,&m) != EOF){
        point[i-1].lng = n;
        point[i-1].lat = m;
    }

	i = 0;
	
	while(fscanf(fp3,"%d",&j) != EOF){
		node[j] = 1;
	}

    
	for(i = 0; i < NODER; i++){
		max_distance = Dijkstra(i,short_node[i],route[i],network);
	}

	for(i = 0; i < NODER; i++){
		if(node[i] == 1){
			while(primary_c.size() != 0){
				primary_c.pop_back();
			}
			x = 0;
			outputfile << "----------------------" << i << "----------------------" << "\n";
			for(j = 0; j < NODER; j++){	
				if(i != j){
					primary_c.push_back(j);
				}

				if(primary_c.size() == 2){
					outputfile << primary_c.front() << "\t" << primary_c.back() << "\n";
					//Dの制約
					min_distance1 = min_geodiversity(i,route,link,point,primary_c.front(),primary_c.back());
					min_distance2 = min_geodiversity(i,route,link,point,primary_c.back(),primary_c.front());
					if(min_distance1 < min_distance2){
						min_distance = min_distance1;
					}else {
						min_distance = min_distance2;
					}

					outputfile << "\t" << min_distance << "\n";
				
					
					if(primary_c.size() == 2){
						primary_c.pop_back();
					}
				}
				
				if(j == NODER - 1){
					x++;
					j = x;
					while(primary_c.size() != 0){
						primary_c.pop_back();
					}
				}
			}
		}
	}	


	
	outputfile.close();

    return 0;
}
