#include "cfd_algo_1.cpp"

int main()
{
    // Generate / Read the boundry Node and Edge information
    // For test we generate a square boundry.
    
    float dimXm=-1,dimYm=-1,dimXM=7,dimYM=3;
    float boundarySpacing = 0.1;

    int NX = static_cast<int>(
        (dimXM - dimXm) / boundarySpacing
    );

    int NY = static_cast<int>(
        (dimYM - dimYm) / boundarySpacing
    );
    int NUM_BOUNDARY_POINTS = 2 * NX + 2 * NY;

    Node* b = new Node[NUM_BOUNDARY_POINTS];

    vector<pair<int,int>> boundaryEdges;
    vector<vector<float>> boundryConstraint;


    // float Uin  = 20.0f;
    float Uin  = 28;
    float Pin  = 0.0f;
    float Pout = 0.0f;

    for(int i = 0; i <= NX; i++)
    {
        int id = i;

        b[id].posn[0] =
            dimXm + i * boundarySpacing;

        b[id].posn[1] =
            dimYM;

        b[id].isBounNode = true;

        // Connect consecutive points
        if(i > 0)
        {
            boundaryEdges.push_back(
                {id - 1, id}
            );
            boundryConstraint.push_back(
                {
                    1,
                    Pin,
                    1,
                    Uin,
                    0.0
                }
            );
        }
    }

    for(int i = 1; i <= NY; i++)
    {
        int id = NX + i;

        b[id].posn[0] =
            dimXM;

        b[id].posn[1] =
            dimYM - i * boundarySpacing;

        b[id].isBounNode = true;

        boundaryEdges.push_back(
            {id - 1, id}
        );

        boundryConstraint.push_back(
            {
                1,
                Pout,
                0,
                0.0f,
                0.0f
            }
        );
    }

    for(int i = 1; i <= NX; i++)
    {
        int id = NX + NY + i;

        b[id].posn[0] =
            dimXM - i * boundarySpacing;

        b[id].posn[1] =
            dimYm;

        b[id].isBounNode = true;

        boundaryEdges.push_back(
            {id - 1, id}
        );

        boundryConstraint.push_back(
            {
                1,
                Pin,
                1,
                Uin,
                0.0f
            }
        );
    }

    for(int i = 1; i < NY; i++)
    {
        int id = NX + NY + NX + i;

        b[id].posn[0] =
            dimXm;

        b[id].posn[1] =
            dimYm + i * boundarySpacing;

        b[id].isBounNode = true;

        boundaryEdges.push_back(
            {id - 1, id}
        );

        boundryConstraint.push_back(
            {
                1,
                Pin,
                1,
                Uin,
                0.0f
            }
        );
    }

    int lastPoint =
        NUM_BOUNDARY_POINTS - 1;

    boundaryEdges.push_back(
        {
            lastPoint,
            0
        }
    );

    boundryConstraint.push_back(
        {
            1,
            Pin,
            1,
            Uin,
            0.0f
        }
    );

    for(int i = 0; i < NUM_BOUNDARY_POINTS; i++)
    {
        b[i].isBounNode = true;
    }

    std::ifstream file("clarky.txt");
    int NUM_INTERIOR_POINT=121;
    float val1, val2;
    Node *in=new Node[NUM_INTERIOR_POINT];
    vector<pair<int,int>> interiorEdges;
    vector<vector<float>> interiorConstraint;
    int ctr=0;
    float A=1e-2;
    float wL=0.4;
    float phi=0.2;
    while (file >> val1 >> val2) {
        in[ctr].posn[0]=val1;
        in[ctr].posn[1]=val2;
        in[ctr].isBounNode=true;
        // in[ctr].posn[1]+=A*sin(2*M_PI*(in[ctr].posn[0]/wL)+phi);
        ctr++;
    }
    file.close();

    for(int i=0;i<NUM_INTERIOR_POINT;i++)
    {
        float theta=(10)*M_PI/180;
        in[i].posn[0]*=3;
        in[i].posn[1]*=3;
        in[i].posn[0]=(cos(theta)*in[i].posn[0]+sin(theta)*in[i].posn[1]);
        in[i].posn[1]=(-sin(theta)*in[i].posn[0]+cos(theta)*in[i].posn[1])+1;
        interiorEdges.push_back({i,(i+1)%NUM_INTERIOR_POINT});
        interiorConstraint.push_back({0,0,1,0,0});
    }
    
    float rC=0.08;
    float rF=0.05;
    Mesh M;
    M.b=b;
    M.in=in;
    M.numBp=NUM_BOUNDARY_POINTS;
    M.numIp=NUM_INTERIOR_POINT;
    M.genPoint(boundaryEdges,boundryConstraint,interiorEdges,interiorConstraint,rC,rF);
    savePoints(M.v,M.NUM_POINT);
    cout<<"Points generated"<<endl;
    M.genTri();
    cout<<M.NUM_CELL<<endl;
    cout<<M.NUM_EDGE<<endl;
    // M.genAdaptMesh(0.06,0.1,2);
    CFD NS;
    NS.m=M;
    NS.NUM_ITER=100;
    NS.dt=1e-2;
    NS.solve_unstructured_NS();
    NS.find_performance();
    cout<<"Lift: "<<NS.Lift<<endl;
    cout<<"Drag: "<<NS.Drag<<endl;
    return 0;
}


    // int NUM_INTERIOR_POINT;
    // float omega=-1.5; // 5 radians per second

    // std::ifstream file("new_points.csv");
    // file>>NUM_INTERIOR_POINT;
    // cout<<NUM_INTERIOR_POINT<<endl;
    // NUM_INTERIOR_POINT+=40;
    // Node *in=new Node[NUM_INTERIOR_POINT];
    // float val1,val2;
    // int ctr=0;
    // char comma;
    // while (file >> val1>>comma>> val2) {
    //     in[ctr].posn[0]=val1;
    //     in[ctr].posn[1]=val2;
    //     in[ctr].isBounNode=true;
    //     ctr++;
    // }

    // file.close();

    // int tmp=ctr;
    // vector<pair<int,int>> interiorEdges;
    // vector<vector<float>> interiorConstraint;
    
    // while(ctr < tmp + 40)
    // {
    //     int local = ctr - tmp;

    //     float theta;
    //     float cx;

    //     if(local < 20)
    //     {
    //         // First circle
    //         theta = 2.0f * M_PI * local / 20.0f;
    //         cx = 0.56f;
    //     }
    //     else
    //     {
    //         // Second circle
    //         int j = local - 20;
    //         theta = 2.0f * M_PI * j / 20.0f;
    //         cx = 2.1f + 0.23f + 0.33f;
    //     }

    //     in[ctr].posn[0] = cx + 0.33f * cos(theta);
    //     in[ctr].posn[1] = 0.3f + 0.33f * sin(theta);

    //     in[ctr].isBounNode = true;

    //     // Connect points within the appropriate circle
    //     if(local < 20)
    //     {
    //         int next = (local + 1) % 20;
    //         interiorEdges.push_back({ctr, tmp + next});
    //     }
    //     else
    //     {
    //         int j = local - 20;
    //         int next = (j + 1) % 20;
    //         interiorEdges.push_back({ctr, tmp + 20 + next});
    //     }

    //     interiorConstraint.push_back(
    //         {0, 0, 1, -Uin * sin(theta), Uin * cos(theta)}
    //     );

    //     ctr++;
    // }


    // for(int i=0;i<NUM_INTERIOR_POINT;i++)
    // {
    //     in[i].posn[0]+=0.3;
    //     in[i].posn[1]-=0.1;
    // }
    // int num_e;
    // std::ifstream file1("new_edges.csv");
    // file1>>num_e;
    // int ae,be;
    // while(file1>>ae>>comma>>be)
    // {
    //     interiorEdges.push_back({ae,be});
    //     interiorConstraint.push_back({0,0,1,0,0});
    // }
    // file1.close();

    // float R=1;
    // float c=3;
    // float thick=0.12;
    // for(int i=0;i<=NUM_INTERIOR_POINT/2;i++)
    // {
    //     if(i==0)
    //     {
    //         in[0].posn[0]=0;
    //         in[0].posn[1]=0;
    //         in[i].isBounNode=true;
    //     }
    //     else if(i==NUM_INTERIOR_POINT/2)
    //     {
    //         in[i].posn[0]=c;
    //         in[i].posn[1]=0;
    //         in[i].isBounNode=true;
    //     }
    //     else
    //     {
    //         float x=c*(1-cos((M_PI*2*i/NUM_INTERIOR_POINT)))/2;
    //         float y=5*thick*c*(0.2969*(sqrt(x/c))-0.126*(x/c)-0.3516*(pow(x/c,2))+0.2843*(pow(x/c,3))-0.1015*(pow(x/c,4)));
    //         in[i].posn[0]=x;
    //         in[i].posn[1]=y;
    //         in[NUM_INTERIOR_POINT-i].posn[0]=x;
    //         in[NUM_INTERIOR_POINT-i].posn[1]=-y;
    //         in[i].isBounNode=true;
    //         in[NUM_INTERIOR_POINT-i].isBounNode=true;
    //     }
    // }
    // for(int i=0;i<NUM_INTERIOR_POINT;i++)
    // {
    //     float theta=(5)*M_PI/180;
    //     in[i].posn[0]=cos(theta)*in[i].posn[0]+sin(theta)*in[i].posn[1];
    //     in[i].posn[1]=-sin(theta)*in[i].posn[0]+cos(theta)*in[i].posn[1];
    //     interiorEdges.push_back({i,(i+1)%NUM_INTERIOR_POINT});
    //     interiorConstraint.push_back({0,0,1,0,0});
    // }
    // for(int i=0;i<NUM_INTERIOR_POINT;i++)
    // {
    //     float theta=2*(float(i)/float(NUM_INTERIOR_POINT))*M_PI;
    //     in[i].posn[0]=R*cos(theta);
    //     in[i].posn[1]=R*sin(theta);
    //     in[i].isBounNode=true;
    //     interiorEdges.push_back({i,(i+1)%NUM_INTERIOR_POINT});
    //     interiorConstraint.push_back({0,0,1,-omega*R*sin(theta),omega*R*cos(theta)});
    // }

    // Node* in = new Node[54];

    // in[0].posn[0]  = 0.000000;
    // in[0].posn[1]  = 0.000000;

    // in[1].posn[0]  = 0.006690;
    // in[1].posn[1]  = 0.011070;

    // in[2].posn[0]  = 0.023980;
    // in[2].posn[1]  = 0.019390;

    // in[3].posn[0]  = 0.073820;
    // in[3].posn[1]  = 0.033720;

    // in[4].posn[0]  = 0.148900;
    // in[4].posn[1]  = 0.046650;

    // in[5].posn[0]  = 0.249170;
    // in[5].posn[1]  = 0.056470;

    // in[6].posn[0]  = 0.349510;
    // in[6].posn[1]  = 0.060300;

    // in[7].posn[0]  = 0.449850;
    // in[7].posn[1]  = 0.058610;

    // in[8].posn[0]  = 0.550130;
    // in[8].posn[1]  = 0.052350;

    // in[9].posn[0]  = 0.650320;
    // in[9].posn[1]  = 0.042640;

    // in[10].posn[0] = 0.750380;
    // in[10].posn[1] = 0.030610;

    // in[11].posn[0] = 0.850300;
    // in[11].posn[1] = 0.017610;

    // in[12].posn[0] = 0.950100;
    // in[12].posn[1] = 0.005300;

    // in[13].posn[0] = 1.000000;
    // in[13].posn[1] = 0.000000;

    // in[14].posn[0] = 0.949900;
    // in[14].posn[1] = 0.001200;

    // in[15].posn[0] = 0.899790;
    // in[15].posn[1] = -0.000870;

    // in[16].posn[0] = 0.799640;
    // in[16].posn[1] = -0.008220;

    // in[17].posn[0] = 0.699640;
    // in[17].posn[1] = -0.017400;

    // in[18].posn[0] = 0.589760;
    // in[18].posn[1] = -0.026440;

    // in[19].posn[0] = 0.500000;
    // in[19].posn[1] = -0.033930;

    // in[20].posn[0] = 0.400320;
    // in[20].posn[1] = -0.038670;

    // in[21].posn[0] = 0.300670;
    // in[21].posn[1] = -0.039660;

    // in[22].posn[0] = 0.200980;
    // in[22].posn[1] = -0.036480;

    // in[23].posn[0] = 0.101180;
    // in[23].posn[1] = -0.028430;

    // in[24].posn[0] = 0.051140;
    // in[24].posn[1] = -0.021210;

    // in[25].posn[0] = 0.013380;
    // in[25].posn[1] = -0.011650;

    // in[26].posn[0] = 0.005700;
    // in[26].posn[1] = -0.007760;

    // for(int j=0;j<27;j++)
    // {
    //     in[j+27]=in[j];
    // }

    // for(int j = 0; j < 27; j++)
    // {
    //     in[j].isBounNode = true;
    //     float theta=(20)*M_PI/180;
    //     float x=cos(theta)*in[j].posn[0]+sin(theta)*in[j].posn[1];
    //     float y=-sin(theta)*in[j].posn[0]+cos(theta)*in[j].posn[1];
    //     in[j].posn[0]=3*x+0.5;
    //     in[j].posn[1]=3*y+0.2;
    // }

    // for(int j = 27; j < 54; j++)
    // {
    //     in[j].isBounNode = true;
    //     float theta=(-10)*M_PI/180;
    //     float x=cos(theta)*in[j].posn[0]+sin(theta)*in[j].posn[1];
    //     float y=-sin(theta)*in[j].posn[0]+cos(theta)*in[j].posn[1];
    //     in[j].posn[0]=2*x+0.1;
    //     in[j].posn[1]=2*y+1.2;
    // }

    // vector<pair<int,int>> interiorEdges ;
    // {
    //     {0,1},
    //     {1,2},
    //     {2,3},
    //     {3,4},
    //     {4,5},
    //     {5,6},
    //     {6,7},
    //     {7,8},
    //     {8,9},
    //     {9,10},
    //     {10,11},
    //     {11,12},
    //     {12,13},

    //     {13,14},
    //     {14,15},
    //     {15,16},
    //     {16,17},
    //     {17,18},
    //     {18,19},
    //     {19,20},
    //     {20,21},
    //     {21,22},
    //     {22,23},
    //     {23,24},
    //     {24,25},
    //     {25,26},

    //     {26,0},

    //     {0+27,1+27},
    //     {1+27,2+27},
    //     {2+27,3+27},
    //     {3+27,4+27},
    //     {4+27,5+27},
    //     {5+27,6+27},
    //     {6+27,7+27},
    //     {7+27,8+27},
    //     {8+27,9+27},
    //     {9+27,10+27},
    //     {10+27,11+27},
    //     {11+27,12+27},
    //     {12+27,13+27},

    //     {13+27,14+27},
    //     {14+27,15+27},
    //     {15+27,16+27},
    //     {16+27,17+27},
    //     {17+27,18+27},
    //     {18+27,19+27},
    //     {19+27,20+27},
    //     {20+27,21+27},
    //     {21+27,22+27},
    //     {22+27,23+27},
    //     {23+27,24+27},
    //     {24+27,25+27},
    //     {25+27,26+27},

    //     {26+27,0+27}
    // };
    // vector<vector<float>> interiorConstraint;

    // for(int i = 0; i < 54; i++)
    // {
    //     // {pDef, P, vDef, Vx, Vy}
    //     interiorConstraint.push_back({0, 0, 1, 0, 0});
    // }