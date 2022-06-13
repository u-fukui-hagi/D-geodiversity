#include <iostream>
#include <fstream>

using namespace std;

#define NODER 50

int main(void){
    int i;
    ofstream outputfile("init_node.txt");
    
    for(i = 0; i < NODER; i++){
        outputfile << i << "\n";
    }

    return 0;
}