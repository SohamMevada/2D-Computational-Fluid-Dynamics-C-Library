#include<bits/stdc++.h>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>
#define TPH_POISSON_IMPLEMENTATION
#include "include/thinks/tph_poisson.h"
using namespace std;

class Node;
class Edge;
class Cell;
class Mesh;

float cross(Node* A, Node* B, Node* P);

class Node{
    public:
    float posn[2];
    bool isBounNode=false;
    bool operator==(Node a);
    bool isInside(vector<Edge*>& e);
    Cell* findCell(list<Cell*>& tri);
};

class Edge{
    public:
    Node* vertex[2]={nullptr,nullptr};
    bool isBounEdge=false;
    bool isInsideEdge=false;
    bool isConstrained=false;
    Cell* neighbor[2]={nullptr,nullptr};
    float vel[2]={1.0f,0.0f};
    bool isVdef=false;
    bool isPdef=false;
    float pDef=0.0f;

    Edge* findEdge(list<Edge*> &edg);
    void findIntersectingEdge(list<Edge*> &edg,queue<Edge*> &intE);
    bool isConvex();
    bool checkIntersect(Edge e);
    bool checkDelaunay();
    void swapEdge(vector<Edge*>& oppEdg,Node *v);
    void findCoeff_momentum_staggered(float *Ax,float *Ay,float *bx,float *by,int row,int NUM_EDGE,float dt,Edge *edg);
};

class Cell{
    public:
    Node* vertex[3]={nullptr,nullptr,nullptr}; //not used till now should we use it or not? we can get the vertices directly from the edges
    Edge* e[3]={nullptr,nullptr,nullptr};
    float P=0.0f;
    float visc=1.81*1e-5;
    float rho=1.24;

    bool inCircumCirc(Node a);
    void findCoeff_continuity(float *A,float *b,int cpos,int NUM_CELL,float dt,Cell *tri);
};

bool Node::operator==(Node a)
{
    bool c1=posn[0]==a.posn[0];
    bool c2=posn[1]==a.posn[1];
    bool c3=isBounNode==a.isBounNode;
    return c1&c2&c3;
}
bool Node::isInside(vector<Edge*>& e)
{

    int i=0;
    while(i<e.size())
    {
        Node *start;
        start=(e[i]->vertex[0]);
        bool inside=false;
        while(i<e.size())
        {
            Node from,to;
            from=*(e[i]->vertex[0]);
            to=*(e[i]->vertex[1]);
            if((from.posn[1]>posn[1])!=(to.posn[1]>posn[1]))
            {  
                float x_inter=from.posn[0]+(posn[1]-from.posn[1])*(to.posn[0]-from.posn[0])/(to.posn[1]-from.posn[1]);
                if(x_inter>posn[0])
                {
                    inside=!inside;
                }
            }
            if(start==e[i]->vertex[1])
            {
                break;
            }
            i++;
        }
        if(inside)
        {
            return true;
        }
        i++;
    }
    return false; 
}
Cell* Node::findCell(list<Cell*>& tri)
{
    for(Cell* c : tri)
    {
        Node* A = c->vertex[0];
        Node* B = c->vertex[1];
        Node* C = c->vertex[2];

        float c1 = cross(A, B, this);
        float c2 = cross(B, C, this);
        float c3 = cross(C, A, this);

        // Point is inside if all have the same sign
        if ((c1 >= 0 && c2 >= 0 && c3 >= 0)||(c1 <= 0 && c2 <= 0 && c3 <= 0))
        {
            return c;
        }
    }
    return nullptr;
}

float cross(Node* A, Node* B, Node* P)
{
    float ABx = B->posn[0] - A->posn[0];
    float ABy = B->posn[1] - A->posn[1];

    float APx = P->posn[0] - A->posn[0];
    float APy = P->posn[1] - A->posn[1];

    return ABx * APy - ABy * APx;
}

Edge* Edge::findEdge(list<Edge*> &edg)
{
    for(auto it=edg.begin();it!=edg.end();++it)
    {
        Edge* e;
        e=*it;
        if((e->vertex[0]==vertex[0] && e->vertex[1]==vertex[1]) || (e->vertex[1]==vertex[0] && e->vertex[0]==vertex[1]))
        {
            return e;
        }
    }
    return nullptr;
}
void Edge::findIntersectingEdge(list<Edge*> &edg,queue<Edge*> &intE)
{
    // populate the queue using the edg list and the edge which is intersecting is the calling class.
    Edge ce;
    ce.vertex[0]=vertex[0];
    ce.vertex[1]=vertex[1];

    for(auto it=edg.begin();it!=edg.end();++it)
    {
        Edge* e;
        e=*it;
        bool intr=ce.checkIntersect(*e);
        if(intr)
        {
            intE.push(e);
        }
    }
}
bool Edge::isConvex()
{
    Cell *c0 = neighbor[0];
    Cell *c1 = neighbor[1];

    Node *A = vertex[0];
    Node *B = vertex[1];

    Node *P = nullptr;
    Node *Q = nullptr;

    // Opposite vertex in c0
    for(int i=0; i<3; i++)
    {
        if(c0->vertex[i] != A && c0->vertex[i] != B)
        {
            P = c0->vertex[i];
            break;
        }
    }

    // Opposite vertex in c1
    for(int i=0; i<3; i++)
    {
        if(c1->vertex[i] != A && c1->vertex[i] != B)
        {
            Q = c1->vertex[i];
            break;
        }
    }

    float sideP = cross(A, B, P);
    float sideQ = cross(A, B, Q);

    float sideA = cross(P, Q, A);
    float sideB = cross(P, Q, B);

    return (sideP * sideQ <= 0) &&
           (sideA * sideB <= 0);
}
bool Edge::checkIntersect(Edge e)
{
    // Check if the 2 edges are equal or not also return false if equal.
    Node* A=vertex[0];
    Node* B=vertex[1];

    Node* C=e.vertex[0];
    Node* D=e.vertex[1];

    float o1=cross(A, B, C);
    float o2=cross(A, B, D);
    float o3=cross(C, D, A);
    float o4=cross(C, D, B);
    // float eps=1e-6;
    return (o1*o2<0.0f && o3*o4<0.0f);
}

bool Edge::checkDelaunay()
{
    Cell *c1 = neighbor[0];
    Cell *c2 = neighbor[1];

    if(c1 == nullptr || c2 == nullptr)
        return true;
    if(!this->isConvex())
        return true;
    Node* n1 = nullptr;

    for(int i = 0; i < 3; i++)
    {
        if(c1->vertex[i] != vertex[0] &&
           c1->vertex[i] != vertex[1])
        {
            n1 = c1->vertex[i];
            break;
        }
    }

    // If opposite point is inside circumcircle,
    // edge is NOT Delaunay.
    bool inside = c2->inCircumCirc(*n1);

    return !inside;
}

void Edge::swapEdge(vector<Edge*>& oppEdg,Node* v)
{
    
    Cell *c1,*c2;
    c1=neighbor[0];
    c2=neighbor[1];

    Node *n1=nullptr,*n2=nullptr,*n3=nullptr,*n4=nullptr;
    n1=vertex[0];
    n2=vertex[1];
    Edge *e11=nullptr,*e21=nullptr,*e12=nullptr,*e22=nullptr;
    bool c1v=false;
    // find n3 and n4
    // also find edges shared by the nodes n3 and n4
    for(int i=0;i<3;i++)
    {
        // for cell 1
        if(!(c1->vertex[i]==n1 || c1->vertex[i]==n2))
        {
            n3=c1->vertex[i];
        }
        if(c1->e[i]!=this)
        {
            if((c1->e[i]->vertex[0]==n1||c1->e[i]->vertex[1]==n1))
            {
                e11=c1->e[i];
            }
            else
            {
                e21=c1->e[i];
            }
        }
        // for cell 2
        if(!(c2->vertex[i]==n1 || c2->vertex[i]==n2))
        {
            n4=c2->vertex[i];
        }
        if(c2->e[i]!=this)
        {
            if((c2->e[i]->vertex[0]==n1||c2->e[i]->vertex[1]==n1))
            {
                e12=c2->e[i];
            }
            else
            {
                e22=c2->e[i];
            }
        }

        // check if v exits in cell 1 or cell 2
        if(c1->vertex[i]==v)
        {
            c1v=true;
        }
    }
    this->vertex[0]=n3;
    this->vertex[1]=n4;

    if(e21->neighbor[0]==c1)
    {
        e21->neighbor[0]=c2;
    }
    else
    {
        e21->neighbor[1]=c2;
    }

    if(e12->neighbor[0]==c2)
    {
        e12->neighbor[0]=c1;
    }
    else
    {
        e12->neighbor[1]=c1;
    }
    
    c1->e[0]=this;
    c1->e[1]=e11;
    c1->e[2]=e12;
    c1->vertex[0]=n1;
    c1->vertex[1]=n3;
    c1->vertex[2]=n4;


    c2->e[0]=this;
    c2->e[1]=e21;
    c2->e[2]=e22;
    c2->vertex[0]=n2;
    c2->vertex[1]=n3;
    c2->vertex[2]=n4;
    
    if(c1v)
    {
        oppEdg.push_back(e12);
        oppEdg.push_back(e22);   
    }
    else
    {
        oppEdg.push_back(e11);
        oppEdg.push_back(e21);
    }

}

bool Cell::inCircumCirc(Node a)
{
    // First Find the Orientation Of Nodes.
    float D_orient;
    D_orient=(vertex[1]->posn[0]-vertex[0]->posn[0])*(vertex[2]->posn[1]-vertex[0]->posn[1])-(vertex[2]->posn[0]-vertex[0]->posn[0])*(vertex[1]->posn[1]-vertex[0]->posn[1]);
    float D_circ;
    float ax=vertex[0]->posn[0]-a.posn[0],ay=vertex[0]->posn[1]-a.posn[1],bx=vertex[1]->posn[0]-a.posn[0],by=vertex[1]->posn[1]-a.posn[1],cx=vertex[2]->posn[0]-a.posn[0],cy=vertex[2]->posn[1]-a.posn[1];
    float a_2=ax*ax+ay*ay,b_2=bx*bx+by*by,c_2=cx*cx+cy*cy;
    D_circ=a_2*(bx*cy-by*cx)-b_2*(ax*cy-ay*cx)+c_2*(ax*by-ay*bx);
    bool inside=false;
    const float EPS=1e-7;
    if(D_orient>EPS)
    {
        // CCW
        inside=D_circ>0;
    }
    else if(D_orient<-EPS)
    {
        // CW
        inside=D_circ<0;
    }

    return inside;
}

void find_normal(float cx,float cy,float ax,float ay,float bx,float by,float *nf)
{
    nf[0]=-(by-ay);
    nf[1]=(bx-ax);
    if((nf[0]*(cx-ax)+nf[1]*(cy-ay))>0)
    {
        nf[0]*=-1;
        nf[1]*=-1;
    }
    float norm=sqrt(nf[0]*nf[0]+nf[1]*nf[1]);
    nf[0]/=norm;
    nf[1]/=norm;
}

void Edge::findCoeff_momentum_staggered(float *Ax,float *Ay,float *bx,float *by,int row,int NUM_EDGE,float dt,Edge *edg)
{
    if(this->isBounEdge)
    {
        Cell *cp;
        cp=(this->neighbor[0]==nullptr)?(this->neighbor[1]):(this->neighbor[0]);
        float xP[2]={0.0f,0.0f};

        for(int i=0;i<3;i++)
        {
            xP[0]+=cp->vertex[i]->posn[0]/3;
            xP[1]+=cp->vertex[i]->posn[1]/3;
        }

        Edge *e_v0_p,*e_v1_p; // Edge sharing vertex 0/1 with cell p
        Node *v0,*v1;
        v0=this->vertex[0];
        v1=this->vertex[1];
        for(int i=0;i<3;i++)
        {
            if((cp->e[i]->vertex[0]==v0 && cp->e[i]->vertex[1]!=v1)||(cp->e[i]->vertex[1]==v0 && cp->e[i]->vertex[0]!=v1))
            {
                e_v0_p=cp->e[i];
            }
            else if((cp->e[i]->vertex[0]!=v0 && cp->e[i]->vertex[1]==v1)||(cp->e[i]->vertex[1]!=v0 && cp->e[i]->vertex[0]==v1))
            {
                e_v1_p=cp->e[i];
            }

        }

        Eigen::Matrix3f mat_p;
        Eigen::Vector3f b_p_x;
        Eigen::Vector3f b_p_y;

        for(int i=0;i<3;i++)
        {
            mat_p(i,0)=1;
            mat_p(i,1)=(cp->e[i]->vertex[0]->posn[0]+cp->e[i]->vertex[1]->posn[0])/2;
            mat_p(i,2)=(cp->e[i]->vertex[0]->posn[1]+cp->e[i]->vertex[1]->posn[1])/2;
            b_p_x[i]=cp->e[i]->vel[0];
            b_p_y[i]=cp->e[i]->vel[1];
        }

        Eigen::Vector3f p_x = mat_p.partialPivLu().solve(b_p_x);
        Eigen::Vector3f p_y = mat_p.partialPivLu().solve(b_p_y);

        float a_p[2],b_p[2],c_p[2];
        a_p[0]=p_x[0]; a_p[1]=p_y[0];
        b_p[0]=p_x[1]; b_p[1]=p_y[1];
        c_p[0]=p_x[2]; c_p[1]=p_y[2];


        float nf1[2],nf2[2],nf3[2];
        find_normal(v1->posn[0],v1->posn[1],v0->posn[0],v0->posn[1],xP[0],xP[1],nf1);
        find_normal(v0->posn[0],v0->posn[1],v1->posn[0],v1->posn[1],xP[0],xP[1],nf2);
        find_normal(xP[0],xP[1],v0->posn[0],v0->posn[1],v1->posn[0],v1->posn[1],nf3);

        float q1[2],q2[2];
        q1[0]=(v0->posn[0]+xP[0])/2;
        q1[1]=(v0->posn[1]+xP[1])/2;
        q2[0]=(v1->posn[0]+xP[0])/2;
        q2[1]=(v1->posn[1]+xP[1])/2;        

        float u_q1[2],u_q2[2];        
        u_q1[0]=a_p[0]+b_p[0]*q1[0]+c_p[0]*q1[1];     
        u_q1[1]=a_p[1]+b_p[1]*q1[0]+c_p[1]*q1[1];     
        u_q2[0]=a_p[0]+b_p[0]*q2[0]+c_p[0]*q2[1];     
        u_q2[1]=a_p[1]+b_p[1]*q2[0]+c_p[1]*q2[1];

        float acell;
        acell=0.5*fabs(v0->posn[0]*(xP[1]-v1->posn[1])+xP[0]*(v1->posn[1]-v0->posn[1])+v1->posn[0]*(v0->posn[1]-xP[1]));

        float flen1,flen2,flen3;
        flen1=sqrt(pow(v0->posn[0]-xP[0],2)+pow(v0->posn[1]-xP[1],2));
        flen2=sqrt(pow(v1->posn[0]-xP[0],2)+pow(v1->posn[1]-xP[1],2));
        flen3=sqrt(pow(v0->posn[0]-v1->posn[0],2)+pow(v0->posn[1]-v1->posn[1],2));
        int epos=this-&edg[0];
        if(this->isVdef)
        {
            // time derivative is 0.
            float flux1=(u_q1[0]*nf1[0]+u_q1[1]*nf1[1])*flen1;
            if(flux1>0)
            {
                Ax[epos+row*NUM_EDGE]+=flux1;
                Ay[epos+row*NUM_EDGE]+=flux1;
            }
            else
            {
                int epos1=e_v0_p-&edg[0];
                Ax[epos1+row*NUM_EDGE]+=flux1;
                Ay[epos1+row*NUM_EDGE]+=flux1;
            }

            float flux2=(u_q2[0]*nf2[0]+u_q2[1]*nf2[1])*flen2;
            if(flux2>0)
            {
                Ax[epos+row*NUM_EDGE]+=flux2;
                Ay[epos+row*NUM_EDGE]+=flux2;
            }
            else
            {
                int epos1=e_v1_p-&edg[0];
                Ax[epos1+row*NUM_EDGE]+=flux2;
                Ay[epos1+row*NUM_EDGE]+=flux2;
            }

            float flux3=(this->vel[0]*nf3[0]+this->vel[1]*nf3[1])*flen3;
            if(flux3<0)
            {
                Ax[epos+row*NUM_EDGE]+=flux3;
                Ay[epos+row*NUM_EDGE]+=flux3;
            }
            else
            {
                Ax[epos+row*NUM_EDGE]+=flux3;
                Ay[epos+row*NUM_EDGE]+=flux3;
            }

            bx[row]+=cp->visc*(b_p[0]*nf1[0]+c_p[0]*nf1[1])*flen1/cp->rho;
            by[row]+=cp->visc*(b_p[1]*nf1[0]+c_p[1]*nf1[1])*flen1/cp->rho;

            bx[row]+=cp->visc*(b_p[0]*nf2[0]+c_p[0]*nf2[1])*flen2/cp->rho;
            by[row]+=cp->visc*(b_p[1]*nf2[0]+c_p[1]*nf2[1])*flen2/cp->rho;

            bx[row]+=cp->visc*(b_p[0]*nf3[0]+c_p[0]*nf3[1])*flen3/cp->rho;
            by[row]+=cp->visc*(b_p[1]*nf3[0]+c_p[1]*nf3[1])*flen3/cp->rho;

        }
        else
        {

            Ax[epos+row*NUM_EDGE]+=(acell/dt);
            Ay[epos+row*NUM_EDGE]+=(acell/dt);

            float flux1=(u_q1[0]*nf1[0]+u_q1[1]*nf1[1])*flen1;
            if(flux1>0)
            {
                Ax[epos+row*NUM_EDGE]+=flux1;
                Ay[epos+row*NUM_EDGE]+=flux1;
            }
            else
            {
                int epos1=e_v0_p-&edg[0];
                Ax[epos1+row*NUM_EDGE]+=flux1;
                Ay[epos1+row*NUM_EDGE]+=flux1;
            }

            float flux2=(u_q2[0]*nf2[0]+u_q2[1]*nf2[1])*flen2;
            if(flux2>0)
            {
                Ax[epos+row*NUM_EDGE]+=flux2;
                Ay[epos+row*NUM_EDGE]+=flux2;
            }
            else
            {
                int epos1=e_v1_p-&edg[0];
                Ax[epos1+row*NUM_EDGE]+=flux2;
                Ay[epos1+row*NUM_EDGE]+=flux2;
            }

            float flux3=(this->vel[0]*nf3[0]+this->vel[1]*nf3[1])*flen3;
            if(flux3<0)
            {
                Ax[epos+row*NUM_EDGE]+=flux3;
                Ay[epos+row*NUM_EDGE]+=flux3;
            }
            else
            {
                Ax[epos+row*NUM_EDGE]+=flux3;
                Ay[epos+row*NUM_EDGE]+=flux3;
            }

            bx[row]+=(acell/dt)*this->vel[0];
            by[row]+=(acell/dt)*this->vel[1];

            bx[row]+=cp->visc*(b_p[0]*nf1[0]+c_p[0]*nf1[1])*flen1/cp->rho;
            by[row]+=cp->visc*(b_p[1]*nf1[0]+c_p[1]*nf1[1])*flen1/cp->rho;

            bx[row]+=cp->visc*(b_p[0]*nf2[0]+c_p[0]*nf2[1])*flen2/cp->rho;
            by[row]+=cp->visc*(b_p[1]*nf2[0]+c_p[1]*nf2[1])*flen2/cp->rho;

            // bx[epos]+=(b_p[0]*nf3[0]+c_p[0]*nf3[1])*flen3;
            // by[epos]+=(b_p[1]*nf3[0]+c_p[1]*nf3[1])*flen3;
        }

    }
    else
    {
        Cell *cn,*cp;
        cp=this->neighbor[0];
        cn=this->neighbor[1];

        float xP[2]={0.0f,0.0f};
        float xN[2]={0.0f,0.0f};

        for(int i=0;i<3;i++)
        {
            xP[0]+=cp->vertex[i]->posn[0]/3;
            xP[1]+=cp->vertex[i]->posn[1]/3;
            xN[0]+=cn->vertex[i]->posn[0]/3;
            xN[1]+=cn->vertex[i]->posn[1]/3;
        }

        Edge *e_v0_p,*e_v1_p,*e_v0_n,*e_v1_n; // Edge sharing vertex 0/1 with cell p/n
        Node *v0,*v1;
        v0=this->vertex[0];
        v1=this->vertex[1];
        for(int i=0;i<3;i++)
        {
            if((cp->e[i]->vertex[0]==v0 && cp->e[i]->vertex[1]!=v1)||(cp->e[i]->vertex[1]==v0 && cp->e[i]->vertex[0]!=v1))
            {
                e_v0_p=cp->e[i];
            }
            else if((cp->e[i]->vertex[0]!=v0 && cp->e[i]->vertex[1]==v1)||(cp->e[i]->vertex[1]!=v0 && cp->e[i]->vertex[0]==v1))
            {
                e_v1_p=cp->e[i];
            }

            if((cn->e[i]->vertex[0]==v0 && cn->e[i]->vertex[1]!=v1)||(cn->e[i]->vertex[1]==v0 && cn->e[i]->vertex[0]!=v1))
            {
                e_v0_n=cn->e[i];
            }
            else if((cn->e[i]->vertex[0]!=v0 && cn->e[i]->vertex[1]==v1)||(cn->e[i]->vertex[1]!=v0 && cn->e[i]->vertex[0]==v1))
            {
                e_v1_n=cn->e[i];
            }
        }

        Eigen::Matrix3f mat_p;
        Eigen::Matrix3f mat_n;
        Eigen::Vector3f b_p_x;
        Eigen::Vector3f b_p_y;
        Eigen::Vector3f b_n_x;
        Eigen::Vector3f b_n_y;

        for(int i=0;i<3;i++)
        {
            mat_p(i,0)=1;
            mat_p(i,1)=(cp->e[i]->vertex[0]->posn[0]+cp->e[i]->vertex[1]->posn[0])/2;
            mat_p(i,2)=(cp->e[i]->vertex[0]->posn[1]+cp->e[i]->vertex[1]->posn[1])/2;
            b_p_x[i]=cp->e[i]->vel[0];
            b_p_y[i]=cp->e[i]->vel[1];
            mat_n(i,0)=1;
            mat_n(i,1)=(cn->e[i]->vertex[0]->posn[0]+cn->e[i]->vertex[1]->posn[0])/2;
            mat_n(i,2)=(cn->e[i]->vertex[0]->posn[1]+cn->e[i]->vertex[1]->posn[1])/2;
            b_n_x[i]=cn->e[i]->vel[0];
            b_n_y[i]=cn->e[i]->vel[1];
        }

        Eigen::Vector3f p_x = mat_p.partialPivLu().solve(b_p_x);
        Eigen::Vector3f p_y = mat_p.partialPivLu().solve(b_p_y);
        Eigen::Vector3f n_x = mat_n.partialPivLu().solve(b_n_x);
        Eigen::Vector3f n_y = mat_n.partialPivLu().solve(b_n_y);

        float a_p[2],b_p[2],c_p[2];
        a_p[0]=p_x[0]; a_p[1]=p_y[0];
        b_p[0]=p_x[1]; b_p[1]=p_y[1];
        c_p[0]=p_x[2]; c_p[1]=p_y[2];
        float a_n[2],b_n[2],c_n[2];
        a_n[0]=n_x[0]; a_n[1]=n_y[0];
        b_n[0]=n_x[1]; b_n[1]=n_y[1];
        c_n[0]=n_x[2]; c_n[1]=n_y[2];

        float nf1[2],nf2[2],nf3[2],nf4[2];
        // face 1 v0-p
        // face 2 v1-p
        // face 3 v0-n
        // face 4 v1-n
        find_normal(v1->posn[0],v1->posn[1],v0->posn[0],v0->posn[1],xP[0],xP[1],nf1);
        find_normal(v0->posn[0],v0->posn[1],v1->posn[0],v1->posn[1],xP[0],xP[1],nf2);
        find_normal(v1->posn[0],v1->posn[1],v0->posn[0],v0->posn[1],xN[0],xN[1],nf3);
        find_normal(v0->posn[0],v0->posn[1],v1->posn[0],v1->posn[1],xN[0],xN[1],nf4);

        float q1[2],q2[2],q3[2],q4[4];
        q1[0]=(v0->posn[0]+xP[0])/2;
        q1[1]=(v0->posn[1]+xP[1])/2;
        q2[0]=(v1->posn[0]+xP[0])/2;
        q2[1]=(v1->posn[1]+xP[1])/2;
        q3[0]=(v0->posn[0]+xN[0])/2;
        q3[1]=(v0->posn[1]+xN[1])/2;
        q4[0]=(v1->posn[0]+xN[0])/2;
        q4[1]=(v1->posn[1]+xN[1])/2;

        float u_q1[2],u_q2[2],u_q3[2],u_q4[2];        
        u_q1[0]=a_p[0]+b_p[0]*q1[0]+c_p[0]*q1[1];     
        u_q1[1]=a_p[1]+b_p[1]*q1[0]+c_p[1]*q1[1];     
        u_q2[0]=a_p[0]+b_p[0]*q2[0]+c_p[0]*q2[1];     
        u_q2[1]=a_p[1]+b_p[1]*q2[0]+c_p[1]*q2[1];

        u_q3[0]=a_n[0]+b_n[0]*q3[0]+c_n[0]*q3[1];     
        u_q3[1]=a_n[1]+b_n[1]*q3[0]+c_n[1]*q3[1];     
        u_q4[0]=a_n[0]+b_n[0]*q4[0]+c_n[0]*q4[1];     
        u_q4[1]=a_n[1]+b_n[1]*q4[0]+c_n[1]*q4[1];

        float acell;
        float sumA,sumB;
        sumA=(v0->posn[0]*xP[1])+(xP[0]*v1->posn[1])+(v1->posn[0]*xN[1])+(xN[0]*v0->posn[1]);
        sumB=(v0->posn[1]*xP[0])+(xP[1]*v1->posn[0])+(v1->posn[1]*xN[0])+(xN[1]*v0->posn[0]);
        acell=0.5*fabs(sumA-sumB);

        float flen1,flen2,flen3,flen4;
        flen1=sqrt(pow(v0->posn[0]-xP[0],2)+pow(v0->posn[1]-xP[1],2));
        flen2=sqrt(pow(v1->posn[0]-xP[0],2)+pow(v1->posn[1]-xP[1],2));
        flen3=sqrt(pow(v0->posn[0]-xN[0],2)+pow(v0->posn[1]-xN[1],2));
        flen4=sqrt(pow(v1->posn[0]-xN[0],2)+pow(v1->posn[1]-xN[1],2));

        int epos=this-&edg[0];

        Ax[epos+row*NUM_EDGE]+=acell/dt;
        Ay[epos+row*NUM_EDGE]+=acell/dt;

        bx[row]+=(acell/dt)*this->vel[0];
        by[row]+=(acell/dt)*this->vel[1];

        float flux1=(u_q1[0]*nf1[0]+u_q1[1]*nf1[1])*flen1;
        if(flux1>0)
        {
            Ax[epos+row*NUM_EDGE]+=flux1;
            Ay[epos+row*NUM_EDGE]+=flux1;
        }
        else
        {
            int epos1=e_v0_p-&edg[0];
            Ax[epos1+row*NUM_EDGE]+=flux1;
            Ay[epos1+row*NUM_EDGE]+=flux1;
        }

        float flux2=(u_q2[0]*nf2[0]+u_q2[1]*nf2[1])*flen2;
        if(flux2>0)
        {
            Ax[epos+row*NUM_EDGE]+=flux2;
            Ay[epos+row*NUM_EDGE]+=flux2;
        }
        else
        {
            int epos1=e_v1_p-&edg[0];
            Ax[epos1+row*NUM_EDGE]+=flux2;
            Ay[epos1+row*NUM_EDGE]+=flux2;
        }

        float flux3=(u_q3[0]*nf3[0]+u_q3[1]*nf3[1])*flen3;
        if(flux3>0)
        {
            Ax[epos+row*NUM_EDGE]+=flux3;
            Ay[epos+row*NUM_EDGE]+=flux3;
        }
        else
        {
            int epos1=e_v0_n-&edg[0];
            Ax[epos1+row*NUM_EDGE]+=flux3;
            Ay[epos1+row*NUM_EDGE]+=flux3;
        }

        float flux4=(u_q4[0]*nf4[0]+u_q4[1]*nf4[1])*flen4;
        if(flux4>0)
        {
            Ax[epos+row*NUM_EDGE]+=flux4;
            Ay[epos+row*NUM_EDGE]+=flux4;
        }
        else
        {
            int epos1=e_v1_n-&edg[0];
            Ax[epos1+row*NUM_EDGE]+=flux4;
            Ay[epos1+row*NUM_EDGE]+=flux4;
        }

        bx[row]+=(cp->visc*flen1/cp->rho)*(b_p[0]*nf1[0]+c_p[0]*nf1[1])+(cp->visc*flen2/cp->rho)*(b_p[0]*nf2[0]+c_p[0]*nf2[1]);
        by[row]+=(cp->visc*flen1/cp->rho)*(b_p[1]*nf1[0]+c_p[1]*nf1[1])+(cp->visc*flen2/cp->rho)*(b_p[1]*nf2[0]+c_p[1]*nf2[1]);

        bx[row]+=(cp->visc*flen3/cp->rho)*(b_n[0]*nf3[0]+c_n[0]*nf3[1])+(cp->visc*flen4/cp->rho)*(b_n[0]*nf4[0]+c_n[0]*nf4[1]);
        by[row]+=(cp->visc*flen3/cp->rho)*(b_n[1]*nf3[0]+c_n[1]*nf3[1])+(cp->visc*flen4/cp->rho)*(b_n[1]*nf4[0]+c_n[1]*nf4[1]);

    }

}

void Cell::findCoeff_continuity(float *A,float *bp,int cpos,int NUM_CELL,float dt,Cell *tri)
{
    float xP[2]={0.0f,0.0f};
    for(int j=0;j<3;j++)
    {
        xP[0]+=(this->e[j]->vertex[0]->posn[0]+this->e[j]->vertex[1]->posn[0])/6;
        xP[1]+=(this->e[j]->vertex[0]->posn[1]+this->e[j]->vertex[1]->posn[1])/6;
    }

    for(int i=0;i<3;i++)
    {
        Edge *e=this->e[i];
        float flen=pow(pow(e->vertex[0]->posn[0]-e->vertex[1]->posn[0],2)+pow(e->vertex[0]->posn[1]-e->vertex[1]->posn[1],2),0.5);
        float fn[2];
        float a[2],b[2],c[2];
        a[0]=e->vertex[0]->posn[0];
        a[1]=e->vertex[0]->posn[1];
        
        b[0]=e->vertex[1]->posn[0];
        b[1]=e->vertex[1]->posn[1];

        int tmp=(this->e[(i+1)%3]->vertex[0]==e->vertex[0] || this->e[(i+1)%3]->vertex[0]==e->vertex[1])?(1):(0);
        c[0]=this->e[(i+1)%3]->vertex[tmp]->posn[0];
        c[1]=this->e[(i+1)%3]->vertex[tmp]->posn[1];
        
        fn[0]=-(b[1]-a[1]);
        fn[1]=b[0]-a[0];
        if((fn[0]*(c[0]-a[0])+fn[1]*(c[1]-a[1]))>0)
        {
            fn[0]*=-1;
            fn[1]*=-1;
        }
        float norm=pow(pow(fn[0],2)+pow(fn[1],2),0.5);
        fn[0]/=norm;
        fn[1]/=norm;
        bp[cpos]+=flen*(this->rho/dt)*(e->vel[0]*fn[0]+e->vel[1]*fn[1]);

        if(e->isPdef&&!e->isVdef)
        {
            //  use the gradient as the Pf to PP.
            float xf[2]={(e->vertex[0]->posn[0]+e->vertex[1]->posn[0])/2,(e->vertex[0]->posn[1]+e->vertex[1]->posn[1])/2};
            float dFP[2]={(xf[0]-xP[0]),(xf[1]-xP[1])};

            // A[cpos+cpos*NUM_CELL]-=(flen/dPf);
            A[cpos+cpos*NUM_CELL]-=flen*((dFP[0]*fn[0]+dFP[1]*fn[1])/(dFP[0]*dFP[0]+dFP[1]*dFP[1]));
            bp[cpos]-=e->pDef*flen*((dFP[0]*fn[0]+dFP[1]*fn[1])/(dFP[0]*dFP[0]+dFP[1]*dFP[1]));//+flen*(e->pDef-this->P)/dPf 
        }
        else if(!e->isBounEdge)
        {
            // use the gradient as the PN to PP.
            Cell *cn=(e->neighbor[0]==this)?(e->neighbor[1]):(e->neighbor[0]);
            float dPN[2];
            float xN[2]={0.0f,0.0f};
            for(int j=0;j<3;j++)
            {
                xN[0]+=(cn->e[j]->vertex[0]->posn[0]+cn->e[j]->vertex[1]->posn[0])/6;
                xN[1]+=(cn->e[j]->vertex[0]->posn[1]+cn->e[j]->vertex[1]->posn[1])/6;
            }

            dPN[0]=xN[0]-xP[0];
            dPN[1]=xN[1]-xP[1];

            int ncpos=cn-&tri[0];
            // A[cpos+cpos*NUM_CELL]-=flen/(dPN[0]*fn[0]+dPN[1]*fn[1]);
            // A[ncpos+cpos*NUM_CELL]+=flen/(dPN[0]*fn[0]+dPN[1]*fn[1]);
            A[cpos+cpos*NUM_CELL]-=flen*(dPN[0]*fn[0]+dPN[1]*fn[1])/(dPN[0]*dPN[0]+dPN[1]*dPN[1]);
            A[ncpos+cpos*NUM_CELL]+=flen*(dPN[0]*fn[0]+dPN[1]*fn[1])/(dPN[0]*dPN[0]+dPN[1]*dPN[1]);
            
            // bp[cpos]+=(flen)*(this->rho/dt)*(e->vel[0]*fn[0]+e->vel[1]*fn[1]); //+(flen)*(cn->P-this->P)/(dPN[0]*fn[0]+dPN[1]*fn[1])

        }
        else
        {
            float xf[2]={(e->vertex[0]->posn[0]+e->vertex[1]->posn[0])/2,(e->vertex[0]->posn[1]+e->vertex[1]->posn[1])/2};
            float dFP[2]={(xf[0]-xP[0]),(xf[1]-xP[1])};

            // A[cpos+cpos*NUM_CELL]-=(flen/dPf);
            // A[cpos+cpos*NUM_CELL]-=flen*((dFP[0]*fn[0]+dFP[1]*fn[1])/(dFP[0]*dFP[0]+dFP[1]*dFP[1]));
            // bp[cpos]+=flen*(this->rho/dt)*(e->vel[0]*fn[0]+e->vel[1]*fn[1]);
        }
    }

}

void dumpEdges1(const Edge *edg, int NUM_EDGE,const char* filename);
void dumpCFDResults(Cell *MC, int NUM_CELL,Edge *ME, int NUM_EDGE,string cname,string ename);

class Mesh{
    public:
    Node* b=nullptr;
    vector<Edge*> bEdge,iEdge;
    int numBp=0;
    Node* in=nullptr;
    int numIp=0;
    vector<Edge*> ce;

    Node* v=nullptr;
    Edge *ME=nullptr;
    Cell *MC=nullptr;
    int NUM_POINT=0;
    int NUM_EDGE=0;
    int NUM_CELL=0;
    
    void genPoint(vector<pair<int,int>>& bedg,vector<vector<float>> &bcns,vector<pair<int,int>>& iedg,vector<vector<float>> &icns,float rc,float rf)
    {
        // Generate Point inside the boundry and outside interior.
        // Our Algorithm For checking Points works as below:
        // Initially we have a list of Nodes and a start node.
        // Then we start performing the check until in the node list we return back to the start node.
        // We continue this process until the list completes.
        // Also we add these points into the list of Nodes V since this is just a input list to the point generator.
        
        Node* t;
        float minX=INFINITY,minY=INFINITY,maxX=-INFINITY,maxY=-INFINITY;
        // t=new Node[numBp+numIp];
        t=new Node[numBp+numIp];
        int num=0;
        // inserting boundry and interior nodes first.
        // inserting boundry points and calculating minimum and maximum coordinates.
        int i=0;
        int ib=0;
        while(ib<numBp)
        {
            t[i]=b[ib];
            if(t[i].posn[0]<minX||t[i].posn[0]>maxX)
            {
                if(t[i].posn[0]<minX)
                    minX=t[i].posn[0];
                else
                    maxX=t[i].posn[0];
            }
            if(t[i].posn[1]<minY||t[i].posn[1]>maxY)
            {
                if(t[i].posn[1]<minY)
                    minY=t[i].posn[1];
                else
                    maxY=t[i].posn[1];
            }
            i++;
            ib++;
        }
        const tph_poisson_real bounds_min[2] ={(tph_poisson_real)(minX+rc),(tph_poisson_real)(minY+rc)};
        const tph_poisson_real bounds_max[2] ={(tph_poisson_real)(maxX-rc),(tph_poisson_real)(maxY-rc)};
        tph_poisson_args argsC = {};
        argsC.bounds_min = bounds_min;
        argsC.bounds_max = bounds_max;
        argsC.radius = (tph_poisson_real)rc;
        argsC.ndims = 2;
        argsC.max_sample_attempts = 30;
        argsC.seed = UINT64_C(1981);
        tph_poisson_sampling samplingC = {};
        int ret_c = tph_poisson_create(&argsC,NULL,&samplingC);
        if(ret_c != TPH_POISSON_SUCCESS)
        {
            std::cout << "Poisson sampling failed. Error: "
                    << ret_c << std::endl;
            return;
        }
        const tph_poisson_real *samplesC =tph_poisson_get_samples(&samplingC);
        if(samplesC == NULL)
        {
            std::cout << "No samples returned." << std::endl;
            tph_poisson_destroy(&samplingC);
            return;
        }


        this->NUM_POINT=numBp+numIp+samplingC.nsamples;

        minX=INFINITY;minY=INFINITY;maxX=-INFINITY;maxY=-INFINITY;
        int ii=0;
        while(ii<numIp)
        {
            t[i]=in[ii];
            if(t[i].posn[0]<minX||t[i].posn[0]>maxX)
            {
                if(t[i].posn[0]<minX)
                    minX=t[i].posn[0];
                else
                    maxX=t[i].posn[0];
            }
            if(t[i].posn[1]<minY||t[i].posn[1]>maxY)
            {
                if(t[i].posn[1]<minY)
                    minY=t[i].posn[1];
                else
                    maxY=t[i].posn[1];
            }
            i++;
            ii++;
        }        

        minX-=3*rf;
        maxX+=3*rf;
        minY-=3*rf;
        maxY+=3*rf;
        const tph_poisson_real bounds_min_i[2] ={(tph_poisson_real)(minX),(tph_poisson_real)(minY)};
        const tph_poisson_real bounds_max_i[2] ={(tph_poisson_real)(maxX),(tph_poisson_real)(maxY)};
        tph_poisson_args argsF = {};
        argsF.bounds_min = bounds_min_i;
        argsF.bounds_max = bounds_max_i;
        argsF.radius = (tph_poisson_real)rf;
        argsF.ndims = 2;
        argsF.max_sample_attempts = 30;
        argsF.seed = UINT64_C(1981);
        tph_poisson_sampling samplingF = {};
        int ret_i = tph_poisson_create(&argsF,NULL,&samplingF);
        if(ret_i != TPH_POISSON_SUCCESS)
        {
            std::cout << "Poisson sampling failed. Error: "
                    << ret_i << std::endl;
            return;
        }
        const tph_poisson_real *samplesF =tph_poisson_get_samples(&samplingF);
        if(samplesF == NULL)
        {
            std::cout << "No samples returned." << std::endl;
            tph_poisson_destroy(&samplingF);
            return;
        }

        this->NUM_POINT+=samplingF.nsamples;

        this->v=new Node[NUM_POINT];
        for(int j=0;j<numBp+numIp;j++)
        {
            v[j]=t[j];
        }
        delete[] t;
        // vector<Edge*> bEdge;
        auto bc=bcns.begin();
        for(auto it=bedg.begin();it!=bedg.end();++it)
        {
            pair<int,int> tmp;
            tmp=*it;
            Edge *e;
            e=new Edge;
            // e->vertex[0]=&t[tmp.first];
            // e->vertex[1]=&t[tmp.second];
            e->vertex[0]=&v[tmp.first];
            e->vertex[1]=&v[tmp.second];
            e->isBounEdge=true;
            e->isInsideEdge=false;

            auto cp=bc->begin();
            e->isPdef=(*(cp++)>0)?(true):(false);
            e->pDef=*(cp++);
            e->isVdef=(*(cp++)>0)?(true):(false);
            e->vel[0]=*(cp++);
            e->vel[1]=*(cp);
            bc++;

            bEdge.push_back(e);
        }
        
        // int ie=0;
        // while(ie<numIp)
        // {
        //     // t[i]=in[ie];
        //     v[i]=in[ie];
        //     i++;
        //     ie++;
        // }

        // vector<Edge*> iEdge;
        auto ic=icns.begin();
        for(auto it=iedg.begin();it!=iedg.end();++it)
        {

            pair<int,int> tmp;
            tmp=*it;
            Edge *e;
            e=new Edge;
            // e->vertex[0]=&t[tmp.first+ib];
            // e->vertex[1]=&t[tmp.second+ib];
            e->vertex[0]=&v[tmp.first+ib];
            e->vertex[1]=&v[tmp.second+ib];
            e->isBounEdge=true;
            e->isInsideEdge=true;

            auto cp=ic->begin();
            e->isPdef=(*(cp++)>0)?(true):(false);
            e->pDef=*(cp++);
            e->isVdef=(*(cp++)>0)?(true):(false);
            e->vel[0]=*(cp++);
            e->vel[1]=*(cp);
            ic++;

            iEdge.push_back(e);
        }

        // concatenate the constraint edges into the ce vector for CDT.
        ce=bEdge;
        ce.reserve(bEdge.size()+iEdge.size());
        ce.insert(ce.end(),iEdge.begin(),iEdge.end());

        // inserting poisson random points into the grid.
        
        for(ptrdiff_t j = 0; j < samplingF.nsamples; j++)
        {
            float x = samplesF[j * samplingF.ndims];
            float y = samplesF[j * samplingF.ndims + 1];
            Node tmp;
            tmp.posn[0]=x;
            tmp.posn[1]=y;
            if(tmp.isInside(bEdge))
            {
                if(!tmp.isInside(iEdge))
                {
                    // t[i]=tmp;
                    v[i]=tmp;
                    i++;
                    if(i==NUM_POINT)
                        break;
                }
            }
        }
        tph_poisson_destroy(&samplingF);

        for(ptrdiff_t j = 0; j < samplingC.nsamples; j++)
        {
            float x = samplesC[j * samplingC.ndims];
            float y = samplesC[j * samplingC.ndims + 1];
            Node tmp;
            tmp.posn[0]=x;
            tmp.posn[1]=y;
            if(!((x>=minX&&x<=maxX)&&(y>=minY&&y<=maxY)))
            {
                if(tmp.isInside(bEdge))
                {
                    v[i]=tmp;
                    i++;
                    if(i==NUM_POINT)
                        break;
                }
            }
            // if(tmp.isInside(bEdge))
            // {
            //     if(!tmp.isInside(iEdge))
            //     {
            //         // t[i]=tmp;
            //         v[i]=tmp;
            //         i++;
            //         if(i==NUM_POINT)
            //             break;
            //     }
            // }
        }
        NUM_POINT=i;
        tph_poisson_destroy(&samplingC);

        // random_device rd;
        // mt19937 g(rd());
        // shuffle(v,v+NUM_POINT,g);
    }

    void genAdaptMesh(float hmin,float hmax,float L)
    {
        int ni=10;
        for(int i=0;i<ni;i++)
        {
            vector<pair<float,float>> force(NUM_POINT,{0.0f,0.0f});
            vector<float> hDist(NUM_POINT,0.0f);

            for(int j=0;j<NUM_POINT;j++)
            {
                float mD=INFINITY;
                for(int k=0;k<numIp;k++)
                {
                    float tmp;
                    tmp=sqrt(pow(v[j].posn[0]-in[k].posn[0],2)+pow(v[j].posn[1]-in[k].posn[1],2));
                    mD=(tmp<mD)?tmp:mD;
                }
                mD=(mD/L>1)?1:mD;
                hDist[j]=hmin+(hmax-hmin)*pow(mD,1);
            }

            float K=50;
            for(int j=0;j<NUM_EDGE;j++)
            {
                int n1,n2;
                Edge *e;
                e=&ME[j];
                n1=e->vertex[0]-&v[0];
                n2=e->vertex[1]-&v[0];
                float meshSize=(hDist[n1]+hDist[n2])/2;
                float edgLen=sqrt(pow(v[n1].posn[0]-v[n2].posn[0],2)+pow(v[n1].posn[1]-v[n2].posn[1],2));
                if(!e->vertex[0]->isBounNode)
                {
                    force[n1].first+=K*(edgLen-meshSize)*(v[n2].posn[0]-v[n1].posn[0]);///edgLen;
                    force[n1].second+=K*(edgLen-meshSize)*(v[n2].posn[1]-v[n1].posn[1]);///edgLen;
                }
                if(!e->vertex[1]->isBounNode)
                {
                    force[n2].first+=-K*(edgLen-meshSize)*(v[n2].posn[0]-v[n1].posn[0]);///edgLen;
                    force[n2].second+=-K*(edgLen-meshSize)*(v[n2].posn[1]-v[n1].posn[1]);///edgLen;
                }
            }

            float stepSize=1e-2;
            vector<pair<float,float>> newPos(NUM_POINT,{0.0f,0.0f});
            for(int j=0;j<NUM_POINT;j++)
            {
                newPos[j].first=v[j].posn[0]+stepSize*force[j].first;
                newPos[j].second=v[j].posn[1]+stepSize*force[j].second;
            }

            for(int j=0;j<NUM_POINT;j++)
            {
                Node tmp;
                tmp.posn[0]=newPos[j].first;
                tmp.posn[1]=newPos[j].second;
                if(tmp.isInside(bEdge))
                {
                    if(!tmp.isInside(iEdge))
                    {
                        v[j].posn[0]=tmp.posn[0];
                        v[j].posn[1]=tmp.posn[1];
                        j++;
                    }
                }
            }

            this->genTri();
        }

    }

    void genTri(){
        // Generate the Edges and Cells using the Constrained Delunay Triangulation.
        // We will first generate the super triangle which is outside the boundry.
        // Then we will do normal delaunay triangulation.
        // Then we will add the constrained edges and replacing intersecting edges.
        // Then we will remove all the cells having the vertices of super triangle.

        // First making the SuperTriangle.
        cout<<"Inside Generate Triangle"<<endl;
        float minX=INFINITY,minY=INFINITY,maxX=-INFINITY,maxY=-INFINITY;
        for(int i=0;i<NUM_POINT;i++)
        {
            if(v[i].posn[0]<minX)
            {
                minX=v[i].posn[0];
            }
            if(v[i].posn[0]>maxX)
            {
                maxX=v[i].posn[0];
            }
            if(v[i].posn[1]<minY)
            {
                minY=v[i].posn[1];
            }
            if(v[i].posn[1]>maxY)
            {
                maxY=v[i].posn[1];
            }
        }
        cout<<minX<<" "<<maxX<<" "<<minY<<" "<<maxY<<endl;

        float a,b;
        minX-=1;
        minY-=1;
        maxX+=1;
        maxY+=1;
        Node* st1=nullptr;
        st1=new Node[3];
        st1[0].posn[0]=(3*minX-maxX)*0.5;
        st1[0].posn[1]=minY;
        st1[1].posn[0]=(3*maxX-minX)*0.5;
        st1[1].posn[1]=minY;
        st1[2].posn[0]=(minX+maxX)*0.5;
        st1[2].posn[1]=2*maxY-minY;

        Edge *se1=nullptr,*se2=nullptr,*se3=nullptr;
        se1=new Edge;
        se2=new Edge;
        se3=new Edge;
        
        se1->vertex[0]=&st1[0];
        se1->vertex[1]=&st1[1];
        se2->vertex[0]=&st1[0];
        se2->vertex[1]=&st1[2];
        se3->vertex[0]=&st1[1];
        se3->vertex[1]=&st1[2];
        
        list<Cell*> tri;
        list<Edge*> edg;
        // We will only play with the above 2 list's in order to generate the CDT.
        edg.push_back(se1);
        edg.push_back(se2);
        edg.push_back(se3);

        Cell* superTri=new Cell;
        superTri->vertex[0]=&st1[0];
        superTri->vertex[1]=&st1[1];
        superTri->vertex[2]=&st1[2];
        superTri->e[0]=se1;
        superTri->e[1]=se2;
        superTri->e[2]=se3;
        tri.push_back(superTri);

        cout<<"Super Triangle Generated"<<endl;
        cout<<st1[0].posn[0]<<" "<<st1[0].posn[1]<<endl;
        cout<<st1[1].posn[0]<<" "<<st1[1].posn[1]<<endl;
        cout<<st1[2].posn[0]<<" "<<st1[2].posn[1]<<endl;

        // Super Triangle generated and stored in list of cells.
        // Now we need to add the new points to generate the delaunay triangulation.
        for(int i=NUM_POINT-1;i>=0;i--)
        // for(int i=0;i<NUM_POINT;i++)
        {
            Cell* st=v[i].findCell(tri);
            // first split the triangle inside which point p lies into 3.
            Edge *e1,*e2,*e3;
            e1=new Edge;
            e2=new Edge;
            e3=new Edge;
            e1->vertex[0]=st->vertex[0];
            e1->vertex[1]=&v[i];
            e2->vertex[0]=st->vertex[1];
            e2->vertex[1]=&v[i];
            e3->vertex[0]=st->vertex[2];
            e3->vertex[1]=&v[i];

            Cell *t1,*t2,*t3;
            t1=new Cell;
            t2=new Cell;
            t3=new Cell;

            e1->neighbor[0]=t1;
            e1->neighbor[1]=t3;
            e2->neighbor[0]=t1;
            e2->neighbor[1]=t2;
            e3->neighbor[0]=t2;
            e3->neighbor[1]=t3;

            t1->vertex[0]=st->vertex[0];
            t1->vertex[1]=st->vertex[1];
            t1->vertex[2]=&v[i];
            t1->e[0]=e1;
            t1->e[1]=e2;
            for(int j=0;j<3;j++)
            {
                bool cond1=st->e[j]->vertex[0]==st->vertex[0] && st->e[j]->vertex[1]==st->vertex[1];
                bool cond2=st->e[j]->vertex[0]==st->vertex[1] && st->e[j]->vertex[1]==st->vertex[0];
                if(cond1 || cond2)
                {
                    t1->e[2]=st->e[j];
                    if(t1->e[2]->neighbor[0]==st)
                    {
                        t1->e[2]->neighbor[0]=t1;
                    }
                    else
                    {
                        t1->e[2]->neighbor[1]=t1;
                    }
                }
            }
            
            t2->vertex[0]=st->vertex[1];
            t2->vertex[1]=st->vertex[2];
            t2->vertex[2]=&v[i];
            t2->e[0]=e2;
            t2->e[1]=e3;
            for(int j=0;j<3;j++)
            {
                bool cond1=st->e[j]->vertex[0]==st->vertex[1] && st->e[j]->vertex[1]==st->vertex[2];
                bool cond2=st->e[j]->vertex[0]==st->vertex[2] && st->e[j]->vertex[1]==st->vertex[1];
                if(cond1 || cond2)
                {
                    t2->e[2]=st->e[j];
                    if(t2->e[2]->neighbor[0]==st)
                    {
                        t2->e[2]->neighbor[0]=t2;
                    }
                    else
                    {
                        t2->e[2]->neighbor[1]=t2;
                    }
                }
            }

            t3->vertex[0]=st->vertex[2];
            t3->vertex[1]=st->vertex[0];
            t3->vertex[2]=&v[i];
            t3->e[0]=e3;
            t3->e[1]=e1;
            for(int j=0;j<3;j++)
            {
                bool cond1=st->e[j]->vertex[0]==st->vertex[0] && st->e[j]->vertex[1]==st->vertex[2];
                bool cond2=st->e[j]->vertex[0]==st->vertex[2] && st->e[j]->vertex[1]==st->vertex[0];
                if(cond1 || cond2)
                {
                    t3->e[2]=st->e[j];
                    if(t3->e[2]->neighbor[0]==st)
                    {
                        t3->e[2]->neighbor[0]=t3;
                    }
                    else
                    {
                        t3->e[2]->neighbor[1]=t3;
                    }
                }
            }
            
            edg.push_back(e1);
            edg.push_back(e2);
            edg.push_back(e3);
            
            tri.push_back(t1);
            tri.push_back(t2);
            tri.push_back(t3);
            
            stack<Edge*> ef;
            ef.push(st->e[0]);
            ef.push(st->e[1]);
            ef.push(st->e[2]);

            tri.remove(st);
            delete st;
            while(!ef.empty())
            {
                Edge* e=ef.top();
                ef.pop();
                Cell *c1,*c2;
                c1=e->neighbor[0];
                c2=e->neighbor[1];
                if(c1==nullptr||c2==nullptr)
                {
                    continue;
                } 
                else
                {
                    bool cD=e->checkDelaunay();
                    if(cD)
                    {
                        continue;
                    }
                    else
                    {
                        // swap the diagonal edge i.e. edge e;
                        vector<Edge*> oppEdg;
                        e->swapEdge(oppEdg,&v[i]);
                        ef.push(oppEdg[0]);
                        ef.push(oppEdg[1]);
                    }
                }
            }
            
        }
        // now we need to add constraint edges. and remove edges which dont follow constraints.
        // suppose we are given the list of constrained edges.

        for(int i=0;i<ce.size();i++)
        {
            Edge* present=ce[i]->findEdge(edg);
            if(present!=nullptr)
            {
                present->isConstrained = true;
                present->isPdef=ce[i]->isPdef;
                present->isVdef=ce[i]->isVdef;
                present->isBounEdge=true;
                present->pDef=ce[i]->pDef;
                present->vel[0]=ce[i]->vel[0];
                present->vel[1]=ce[i]->vel[1];
                present->isInsideEdge=ce[i]->isInsideEdge;
                continue;
            }
            else
            {
                vector<Edge*> crtdEdg;
                bool F=false;
                bool EdgFlp=false;
                while(!F)
                {
                    EdgFlp=false;
                    for(auto it=edg.begin();it!=edg.end();++it)
                    {
                        Edge* e;
                        e=*it;
                        bool intersect=false;
                        intersect=e->checkIntersect(*ce[i]);
                        if(intersect)
                        {
                            if(e->isConvex())
                            {
                                EdgFlp=true;
                                vector<Edge*> dustbin;
                                e->swapEdge(dustbin, e->vertex[0]);
                                crtdEdg.push_back(e);
                                F=(e->vertex[0]==ce[i]->vertex[0] && e->vertex[1]==ce[i]->vertex[1]) || (e->vertex[1]==ce[i]->vertex[0] && e->vertex[0]==ce[i]->vertex[1]);
                                if(F)
                                {
                                    e->isConstrained=true;
                                    e->isPdef=ce[i]->isPdef;
                                    e->isVdef=ce[i]->isVdef;
                                    e->isBounEdge=true;
                                    e->pDef=ce[i]->pDef;
                                    e->vel[0]=ce[i]->vel[0];
                                    e->vel[1]=ce[i]->vel[1];
                                    e->isInsideEdge=ce[i]->isInsideEdge;
                                }
                            }
                        }
                    }
                    if(!EdgFlp)
                    {
                        break;
                    }
                }

                bool swaps=true;
                // We loop over created edges and then 
                // while(swaps)
                // {
                //     swaps=false;
                //     for(int j=0;j<crtdEdg.size();j++)
                //     {
                //         Edge* e1;
                //         e1=crtdEdg[j];
                //         Cell *c1,*c2;
                //         c1=e1->neighbor[0];
                //         c2=e1->neighbor[1];
                //         bool condn=e1->checkDelaunay();
                //         if(!condn&&!e1->isConstrained) //
                //         {
                //             // swap the edge
                //             vector<Edge*> dustbin;
                //             e1->swapEdge(dustbin,e1->vertex[0]);
                //             // crtdEdg.erase(crtdEdg.begin()+j);
                //             // j--;
                //             swaps=true;
                //         }
                //     }
                // }
            }

        }
        // Now we remove the super triangle nodes and cells.
        for(auto it=edg.begin();it!=edg.end();)
        {
            // check if an edge is not in constrained edge but has both vertices as boundry nodes.
            Edge* e=*it;
            Node *a,*b;
            a=e->vertex[0];
            b=e->vertex[1];
            // first find if both vertices belong to Node *in
            if(a->isBounNode && b->isBounNode)
            {
                bool inA=false,inB=false;
                for(int i=0;i<numIp;i++)
                {
                    if(in[i]==*a)
                    {
                        inA=true;
                    }
                    if(in[i]==*b)
                    {
                        inB=true;
                    }
                }
                if(inA&&inB)
                {
                    // now we check if the edge is not present in ce
                    bool isThere=false;
                    for(auto it1=ce.begin();it1!=ce.end();++it1)
                    {
                        Edge *tmp;
                        tmp=*it1;
                        if((tmp->vertex[0]==a && tmp->vertex[1]==b)||(tmp->vertex[1]==a && tmp->vertex[0]==b))
                        {
                            isThere=true;
                            break;
                        }
                    }
                    if(!isThere)
                    {
                        it=edg.erase(it);
                        // clean up
                        Cell *t1,*t2;
                        t1=e->neighbor[0];
                        t2=e->neighbor[1];
                        if(t1!=nullptr)
                        {
                            for(int j=0;j<3;j++)
                            {
                                if(t1->e[j]->neighbor[0]==t1)
                                {
                                    t1->e[j]->neighbor[0]=nullptr;
                                }
                                else if(t1->e[j]->neighbor[1]==t1)
                                {
                                    t1->e[j]->neighbor[1]=nullptr;
                                }
                            }
                            tri.remove(t1);
                            delete t1;
                        }
                        if(t2!=nullptr)
                        {
                            for(int j=0;j<3;j++)
                            {
                                if(t2->e[j]->neighbor[0]==t2)
                                {
                                    t2->e[j]->neighbor[0]=nullptr;
                                }
                                else if(t2->e[j]->neighbor[1]==t2)
                                {
                                    t2->e[j]->neighbor[1]=nullptr;
                                }
                            }
                            tri.remove(t2);
                            delete t2;
                        }
                        
                        delete e;
                    }
                    else
                    {
                        ++it;
                    }
                }
                else
                {
                    ++it;
                }
            }
            else
            {
                ++it;
            }
        }

        for(auto it=tri.begin();it!=tri.end();)
        {
            Cell* c;
            c=*it;
            // check if the cell is a part of super triangle
            bool superC=false;
            for(int j=0;j<3;j++)
            {
                Node* nc;
                nc=c->vertex[j];
                if(nc==&st1[0]||nc==&st1[1]||nc==&st1[2])
                {
                    superC=true;
                }
            }

            if(superC)
            {
                for(int j=0;j<3;j++)
                {
                    if(c->e[j]->neighbor[0]==c)
                    {
                        c->e[j]->neighbor[0]=nullptr;
                    }
                    else if(c->e[j]->neighbor[1]==c)
                    {
                        c->e[j]->neighbor[1]=nullptr;
                    }
                    else{
                        cout<<"Report"<<endl;
                    }
                }
                it=tri.erase(it);
                delete c;
            }
            else
            {
                ++it;
            }

        }
        // now we need to delete the remaining edges
        for(auto it=edg.begin(); it!=edg.end();)
        {
            Edge* ne=*it;

            if(ne->neighbor[0]==nullptr && ne->neighbor[1]==nullptr)
            {
                it=edg.erase(it);
                delete ne;
            }
            else
            {
                ++it;
            }
        }
        // dumpEdges(edg,"edges_after_cleanup.txt");

        // Now we have a constrained Delaunay Triangulated Mesh.
        NUM_EDGE=edg.size();
        NUM_CELL=tri.size();
        // Edge *ME;
        ME=new Edge[NUM_EDGE];
        // Cell *MC;
        MC=new Cell[NUM_CELL];
        map<Edge*,Edge*> edgeLUT;
        map<Cell*,Cell*> cellLUT;
        int j=0;

        for(auto it=edg.begin();it!=edg.end();++it)
        {
            Edge* e;
            e=*it;
            ME[j]=*e;
            edgeLUT[e]=&ME[j];
            j++;
        }
        j=0;
        for(auto it=tri.begin();it!=tri.end();++it)
        {
            Cell* c;
            c=*it;
            MC[j]=*c;
            cellLUT[c]=&MC[j];
            j++;
        }
        for(int i=0;i<NUM_EDGE;i++)
        {
            if(ME[i].neighbor[0]!=nullptr)
                ME[i].neighbor[0]=cellLUT[ME[i].neighbor[0]];
            if(ME[i].neighbor[1]!=nullptr)
                ME[i].neighbor[1]=cellLUT[ME[i].neighbor[1]];
        }
        for(int i=0;i<NUM_CELL;i++)
        {
            MC[i].e[0]=edgeLUT[MC[i].e[0]];
            MC[i].e[1]=edgeLUT[MC[i].e[1]];
            MC[i].e[2]=edgeLUT[MC[i].e[2]];
        }
        dumpEdges1(ME,NUM_EDGE,"edges_after_cleanup.txt");
    }
};

void solve_momentum(float *Ax,float *Ay,float *bx,float *by,Edge *edg,int nr,int ne)
{
    // here we solve the constrained linear equation Au=b using BLAS or CUDA library using GPU/(All CPU cores).
    // The result of solve will be inserted into edg 
    // Firstly we seperate the defined and non defined edges.
    bool add[ne]={false};
    int num_free=0;
    vector<int> MapBack;
    for(int i=0;i<ne;i++)
    {
        if(edg[i].isVdef)
        {
            add[i]=false;
        }
        else
        {
            add[i]=true;
            MapBack.push_back(i);
            num_free++;
        }
    }

    for(int i=0;i<nr*ne;i++)
    {
        int pc=int(i/ne);
        int pe=int(i%ne);
        if(!add[pe])
        {
            bx[pc]-=Ax[i]*edg[pe].vel[0];
            by[pc]-=Ay[i]*edg[pe].vel[1];
        }
    }

    Eigen::SparseMatrix<double> A_x(nr,num_free);
    Eigen::SparseMatrix<double> A_y(nr,num_free);
    Eigen::VectorXd b_x(nr);
    Eigen::VectorXd b_y(nr);
    cout<<nr<<" "<<num_free<<endl;
    for(int r=0;r<nr;r++)
    {
        for(int c=0;c<num_free;c++)
        {
            if(Ax[r*ne+MapBack[c]]!=0)
                A_x.insert(r,c)=Ax[r*ne+MapBack[c]];
            if(Ay[r*ne+MapBack[c]]!=0)
                A_y.insert(r,c)=Ay[r*ne+MapBack[c]];
        }
    }
    for(int i=0;i<nr;i++)
    {
        b_x(i)=bx[i];
        b_y(i)=by[i];
    }

    // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solve_x;
    // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solve_y;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solve_x;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solve_y;

    solve_x.compute(A_x);
    solve_y.compute(A_y);

    Eigen::VectorXd ux = solve_x.solve(b_x);
    Eigen::VectorXd uy = solve_y.solve(b_y);

    cout << "Momentum X solve: "
     << solve_x.info() << endl;

    cout << "Momentum Y solve: "
     << solve_y.info() << endl;

    float maxXvel=0;
    float maxYvel=0;
    float maxUpdtNorm=0;
    for(int i=0;i<num_free;i++)
    {
        if(maxXvel<fabs(edg[MapBack[i]].vel[0]-ux[i]))
            maxXvel=fabs(edg[MapBack[i]].vel[0]-ux[i]);
        if(maxYvel<fabs(edg[MapBack[i]].vel[1]-uy[i]))
            maxYvel=fabs(edg[MapBack[i]].vel[1]-uy[i]);
        float updt=sqrt(pow(edg[MapBack[i]].vel[0]-ux[i],2)+pow(edg[MapBack[i]].vel[1]-uy[i],2));
        if(updt>maxUpdtNorm)
            maxUpdtNorm=updt;
        edg[MapBack[i]].vel[0]=ux[i];
        edg[MapBack[i]].vel[1]=uy[i];
    }
    cout<<"Convergence?: "<<maxUpdtNorm<<endl;
    cout<<"Velocity update x : "<<maxXvel<<endl;
    cout<<"Velocity update y : "<<maxYvel<<endl;
    // Then we use Eigen to firstly convert our dense matrix into a sparse matrix in Eigen
    // Then we solve that sparse matrix Equation using iterative solvers since direct inversion not possible.
    // Then we convert the result back into the *edg format and update its values.
}

void solve_continuity(float *A,float *b,Cell *tri,Edge *edg,int nc,int ne,float dt)
{
    Eigen::SparseMatrix<double> A_p(nc,nc);
    Eigen::VectorXd b_p(nc);
    for(int i=0;i<nc*nc;i++)
    {
        int r=int(i/nc);
        int c=int(i%nc);
        if(A[i]!=0)
        {
            A_p.insert(r,c)=A[i];
            // A_p.insert(r,c,A[i]);
        }
    }
    for(int i=0;i<nc;i++)
    {
        b_p[i]=b[i];
    }
    
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,Eigen::IncompleteLUT<double>> solver_p;
    // solver_p.compute(A_p);
    // Eigen::VectorXd P_n1=solver_p.solve(b_p);
    Eigen::BiCGSTAB<
        Eigen::SparseMatrix<double>,
        Eigen::IncompleteLUT<double>
    > solver_p;

    solver_p.setTolerance(1e-10);
    solver_p.setMaxIterations(5000);

    solver_p.compute(A_p);
    Eigen::VectorXd P_n1 = solver_p.solve(b_p);

    // Velocity correction using continuity equation
    for(int i=0;i<ne;i++)
    {
        Edge *e=&edg[i];
        if(!e->isVdef)
        {
            float fn[2];
            // to find the face normal: (a,b) --> (b,-a) || (-b,a)
            // from these 2 options we need to find the one which is outward facing compared to the cell.
            Cell *cn=(e->neighbor[0]==nullptr)?(e->neighbor[1]):(e->neighbor[0]);

            float a[2],b[2],c[2];
            a[0]=e->vertex[0]->posn[0];
            a[1]=e->vertex[0]->posn[1];
            
            b[0]=e->vertex[1]->posn[0];
            b[1]=e->vertex[1]->posn[1];
            float flen=pow(pow(a[0]-b[0],2)+pow(a[1]-b[1],2),0.5); // added flen since in the FVM momentum equation it is there.
            int k=0;
            for(;k<3;k++)
            {
                if(cn->e[k]==e)
                    break;
            }
            int tmp=(cn->e[(k+1)%3]->vertex[0]==e->vertex[0] || cn->e[(k+1)%3]->vertex[0]==e->vertex[1])?(1):(0);
            c[0]=cn->e[(k+1)%3]->vertex[tmp]->posn[0];
            c[1]=cn->e[(k+1)%3]->vertex[tmp]->posn[1];
            
            fn[0]=-(b[1]-a[1]);
            fn[1]=b[0]-a[0];
            if((fn[0]*(c[0]-a[0])+fn[1]*(c[1]-a[1]))>0)
            {
                fn[0]*=-1;
                fn[1]*=-1;
            }
            float norm=pow(pow(fn[0],2)+pow(fn[1],2),0.5);
            fn[0]/=norm;
            fn[1]/=norm;

            if(e->isBounEdge)
            {
                // here we take the gradient of the face pressure and cell pressure
                Cell *c1;
                c1=(e->neighbor[0]==nullptr)?(e->neighbor[1]):(e->neighbor[0]);
                float dFP[2]={0.0f,0.0f};
                // dFP f --> P
                for(int j=0;j<3;j++)
                {
                    dFP[0]+=(c1->e[j]->vertex[0]->posn[0]+c1->e[j]->vertex[1]->posn[0])/6;
                    dFP[1]+=(c1->e[j]->vertex[0]->posn[1]+c1->e[j]->vertex[1]->posn[1])/6;
                }
                dFP[0]=(e->vertex[0]->posn[0]+e->vertex[1]->posn[0])/2-dFP[0];
                dFP[1]=(e->vertex[0]->posn[1]+e->vertex[1]->posn[1])/2-dFP[1];
                int pc1=c1-&tri[0];
                if(e->isPdef)
                {
                    e->vel[0]-=(dt/c1->rho)*((e->pDef-P_n1[pc1])/(dFP[0]*dFP[0]+dFP[1]*dFP[1]))*dFP[0];
                    e->vel[1]-=(dt/c1->rho)*((e->pDef-P_n1[pc1])/(dFP[0]*dFP[0]+dFP[1]*dFP[1]))*dFP[1];
                }
                // else
                // {
                //     e->vel[0]-=(dt/c1->rho)*((-c1->P)/(dFP[0]*dFP[0]+dFP[1]*dFP[1]))*dFP[0];
                //     e->vel[1]-=(dt/c1->rho)*((-c1->P)/(dFP[0]*dFP[0]+dFP[1]*dFP[1]))*dFP[1];
                // }
            }
            else
            {
                // here we take the gradient of the 2 neighbouring cell pressure.
                Cell *c1,*c2;
                c1=cn;
                c2=(e->neighbor[0]==cn)?(e->neighbor[1]):(e->neighbor[0]);
                int pc1=c1-&tri[0];
                int pc2=c2-&tri[0];
                float d12[2]={0.0f,0.0f};
                // d12 c1 --> c2
                for(int j=0;j<3;j++)
                {
                    d12[0]+=(c2->e[j]->vertex[0]->posn[0]+c2->e[j]->vertex[1]->posn[0])/6-(c1->e[j]->vertex[0]->posn[0]+c1->e[j]->vertex[1]->posn[0])/6;
                    d12[1]+=(c2->e[j]->vertex[0]->posn[1]+c2->e[j]->vertex[1]->posn[1])/6-(c1->e[j]->vertex[0]->posn[1]+c1->e[j]->vertex[1]->posn[1])/6;
                }

                e->vel[0]-=(dt/c1->rho)*((P_n1[pc2]-P_n1[pc1])/(d12[0]*d12[0]+d12[1]*d12[1]))*d12[0];
                e->vel[1]-=(dt/c1->rho)*((P_n1[pc2]-P_n1[pc1])/(d12[0]*d12[0]+d12[1]*d12[1]))*d12[1];

            }
        }
    }


    float alpha=1; // keep 1 to check the continuity now working.
    for(int i=0;i<nc;i++)
    {
        // cout<<tri[i].P<<" "<<P_n1[i]<<endl;
        tri[i].P=(1-alpha)*tri[i].P+alpha*P_n1[i];
    }
}

class CFD{
    public:
    Mesh m;
    float Lift,Drag;
    // float rho=1e2;
    int NUM_ITER=0;
    float dt=1e-1;
    void solve_unstructured_NS()
    {
        for(int i=0;i<NUM_ITER;i++)
        {
            cout<<"Iteration "<<i<<endl;
            solve_NS();
            string cname;
            string ename;
            ostringstream cs;
            ostringstream es;
            cs << "results/cell_data_" << setfill('0') << setw(4) << i << ".csv";
            es << "results/edge_data_" << setfill('0') << setw(4) << i << ".csv";
            cname = cs.str();
            ename = es.str();
            // cname="/results/"+cname;
            // ename="/results/"+ename;
            
            dumpCFDResults(m.MC,m.NUM_CELL,m.ME,m.NUM_EDGE,cname,ename);
        }
        dumpCFDResults(m.MC,m.NUM_CELL,m.ME,m.NUM_EDGE,"cell_data.csv","edge_data.csv");
    }
    void solve_NS()
    {
        // first estimate velocity using momentum
        // build the Ax and Ay matrices also the bx and by
        // to make these equations easy to solve we can A*edgelist=b
        float *Ax=new float[m.NUM_EDGE*m.NUM_EDGE];
        float *Ay=new float[m.NUM_EDGE*m.NUM_EDGE];
        float *bx=new float[m.NUM_EDGE];
        float *by=new float[m.NUM_EDGE];
        for(int i=0;i<m.NUM_EDGE*m.NUM_EDGE;i++)
        {
            Ax[i]=0;
            Ay[i]=0;
        }
        for(int i=0;i<m.NUM_EDGE;i++)
        {
            bx[i]=0;
            by[i]=0;
        }

        int row=0;
        for(int i=0;i<m.NUM_EDGE;i++)
        {
            Edge *e;
            e=&m.ME[i];
            // if(e->vertex[0]->posn[0]==e->vertex[1]->posn[0] && e->vertex[0]->posn[1]==e->vertex[1]->posn[1])
            // {
            //     cout<<"Found a degenerate edge"<<endl;
            //     cout<<e->vertex[0]<<" "<<e->vertex[1]<<endl;
            // }
            if(!e->isVdef){
                // if(row==25 ||row==27||row==998||row==999||row==1000)
                // {
                //     cout<<"Edge posn: "<<e->vertex[0]->posn[0]<<" "<<e->vertex[0]->posn[1]<<" "<<e->vertex[1]->posn[0]<<" "<<e->vertex[1]->posn[1]<<endl;
                // }
                e->findCoeff_momentum_staggered(Ax,Ay,bx,by,row,m.NUM_EDGE,dt,m.ME);
                row++;
            }
        }
        // saveMatrix("a1.csv",Ax,row,m.NUM_EDGE);
        // saveMatrix("b1.csv",bx,row,1);
        solve_momentum(Ax,Ay,bx,by,m.ME,row,m.NUM_EDGE);

        delete[] Ax;
        delete[] Ay;
        delete[] bx;
        delete[] by;
        // then estimate pressure using conservation and update velocity
        float *A=new float[m.NUM_CELL*m.NUM_CELL];
        float *b=new float[m.NUM_CELL];
        for(int i=0;i<m.NUM_CELL*m.NUM_CELL;i++)
        {
            A[i]=0;
        }
        for(int i=0;i<m.NUM_CELL;i++)
        {
            b[i]=0;
        }
        
        for(int i=0;i<m.NUM_CELL;i++)
        {
            Cell *c;
            c=&m.MC[i];
            c->findCoeff_continuity(A,b,i,m.NUM_CELL,dt,m.MC);
        }

        float maxFlux=0;
        float sumFlux=0;
        for(int i=0;i<m.NUM_CELL;i++)
        {
            float flux=0.0f;
            Cell *cp=&m.MC[i];
            for(int j=0;j<3;j++)
            {
                Edge *e=cp->e[j];
                float fn[2];
                float a[2],b[2],c[2];
                a[0]=e->vertex[0]->posn[0];
                a[1]=e->vertex[0]->posn[1];
                
                b[0]=e->vertex[1]->posn[0];
                b[1]=e->vertex[1]->posn[1];
                float flen=pow(pow(a[0]-b[0],2)+pow(a[1]-b[1],2),0.5); // added flen since in the FVM momentum equation it is there.
                
                int tmp=(cp->e[(j+1)%3]->vertex[0]==e->vertex[0] || cp->e[(j+1)%3]->vertex[0]==e->vertex[1])?(1):(0);
                c[0]=cp->e[(j+1)%3]->vertex[tmp]->posn[0];
                c[1]=cp->e[(j+1)%3]->vertex[tmp]->posn[1];
                
                fn[0]=-(b[1]-a[1]);
                fn[1]=b[0]-a[0];
                if((fn[0]*(c[0]-a[0])+fn[1]*(c[1]-a[1]))>0)
                {
                    fn[0]*=-1;
                    fn[1]*=-1;
                }
                float norm=pow(pow(fn[0],2)+pow(fn[1],2),0.5);
                fn[0]/=norm;
                fn[1]/=norm;

                flux+=(e->vel[0]*fn[0]+e->vel[1]*fn[1])*flen;
            }
            sumFlux+=flux;
            // if((cp->e[0]->isBounEdge||cp->e[1]->isBounEdge||cp->e[2]->isBounEdge))
            // {
            //     cout<<"Outside Cell flux: "<<flux<<endl;
            // }
            if(maxFlux<fabs(flux))
            {
                maxFlux=fabs(flux);
            }
        }
        cout<<maxFlux<<endl;
        cout<<sumFlux<<endl;

        solve_continuity(A,b,m.MC,m.ME,m.NUM_CELL,m.NUM_EDGE,dt);
        delete[] A;
        delete[] b;

        maxFlux=0;
        sumFlux=0;
        for(int i=0;i<m.NUM_CELL;i++)
        {
            float flux=0.0f;
            Cell *cp=&m.MC[i];
            for(int j=0;j<3;j++)
            {
                Edge *e=cp->e[j];
                float fn[2];
                float a[2],b[2],c[2];
                a[0]=e->vertex[0]->posn[0];
                a[1]=e->vertex[0]->posn[1];
                
                b[0]=e->vertex[1]->posn[0];
                b[1]=e->vertex[1]->posn[1];
                float flen=pow(pow(a[0]-b[0],2)+pow(a[1]-b[1],2),0.5); // added flen since in the FVM momentum equation it is there.
                
                int tmp=(cp->e[(j+1)%3]->vertex[0]==e->vertex[0] || cp->e[(j+1)%3]->vertex[0]==e->vertex[1])?(1):(0);
                c[0]=cp->e[(j+1)%3]->vertex[tmp]->posn[0];
                c[1]=cp->e[(j+1)%3]->vertex[tmp]->posn[1];
                
                fn[0]=-(b[1]-a[1]);
                fn[1]=b[0]-a[0];
                if((fn[0]*(c[0]-a[0])+fn[1]*(c[1]-a[1]))>0)
                {
                    fn[0]*=-1;
                    fn[1]*=-1;
                }
                float norm=pow(pow(fn[0],2)+pow(fn[1],2),0.5);
                fn[0]/=norm;
                fn[1]/=norm;

                flux+=(e->vel[0]*fn[0]+e->vel[1]*fn[1])*flen;
            }
            sumFlux+=flux;
            // if((cp->e[0]->isBounEdge||cp->e[1]->isBounEdge||cp->e[2]->isBounEdge))
            // {
            //     cout<<"Outside Cell flux: "<<flux<<endl;
            // }
            if(maxFlux<fabs(flux))
            {
                maxFlux=fabs(flux);
            }
        }
        cout<<maxFlux<<endl;
        cout<<sumFlux<<endl;
        // Velocity correction using continuity equation
        // for(int i=0;i<m.NUM_EDGE;i++)
        // {
        //     Edge *e=&m.ME[i];
        //     if(!e->isVdef)
        //     {
        //         float fn[2];
        //         // to find the face normal: (a,b) --> (b,-a) || (-b,a)
        //         // from these 2 options we need to find the one which is outward facing compared to the cell.
        //         Cell *cn=(e->neighbor[0]==nullptr)?(e->neighbor[1]):(e->neighbor[0]);

        //         float a[2],b[2],c[2];
        //         a[0]=e->vertex[0]->posn[0];
        //         a[1]=e->vertex[0]->posn[1];
                
        //         b[0]=e->vertex[1]->posn[0];
        //         b[1]=e->vertex[1]->posn[1];
        //         float flen=pow(pow(a[0]-b[0],2)+pow(a[1]-b[1],2),0.5); // added flen since in the FVM momentum equation it is there.
        //         int k=0;
        //         for(;k<3;k++)
        //         {
        //             if(cn->e[k]==e)
        //                 break;
        //         }
        //         int tmp=(cn->e[(k+1)%3]->vertex[0]==e->vertex[0] || cn->e[(k+1)%3]->vertex[0]==e->vertex[1])?(1):(0);
        //         c[0]=cn->e[(k+1)%3]->vertex[tmp]->posn[0];
        //         c[1]=cn->e[(k+1)%3]->vertex[tmp]->posn[1];
                
        //         fn[0]=-(b[1]-a[1]);
        //         fn[1]=b[0]-a[0];
        //         if((fn[0]*(c[0]-a[0])+fn[1]*(c[1]-a[1]))>0)
        //         {
        //             fn[0]*=-1;
        //             fn[1]*=-1;
        //         }
        //         float norm=pow(pow(fn[0],2)+pow(fn[1],2),0.5);
        //         fn[0]/=norm;
        //         fn[1]/=norm;

        //         if(e->isBounEdge)
        //         {
        //             // here we take the gradient of the face pressure and cell pressure
        //             Cell *c1;
        //             c1=(e->neighbor[0]==nullptr)?(e->neighbor[1]):(e->neighbor[0]);
        //             float dFP[2]={0.0f,0.0f};
        //             // dFP f --> P
        //             for(int j=0;j<3;j++)
        //             {
        //                 dFP[0]+=(c1->e[j]->vertex[0]->posn[0]+c1->e[j]->vertex[1]->posn[0])/6;
        //                 dFP[1]+=(c1->e[j]->vertex[0]->posn[1]+c1->e[j]->vertex[1]->posn[1])/6;
        //             }
        //             dFP[0]=(e->vertex[0]->posn[0]+e->vertex[1]->posn[0])/2-dFP[0];
        //             dFP[1]=(e->vertex[0]->posn[1]+e->vertex[1]->posn[1])/2-dFP[1];
        //             if(e->isPdef)
        //             {
        //                 e->vel[0]-=(dt/rho)*((e->pDef-c1->P)/(dFP[0]*dFP[0]+dFP[1]*dFP[1]))*dFP[0];
        //                 e->vel[1]-=(dt/rho)*((e->pDef-c1->P)/(dFP[0]*dFP[0]+dFP[1]*dFP[1]))*dFP[1];
        //             }
        //             // else
        //             // {
        //             //     e->vel[0]-=(dt/rho)*((-c1->P)/(dFP[0]*dFP[0]+dFP[1]*dFP[1]))*dFP[0];
        //             //     e->vel[1]-=(dt/rho)*((-c1->P)/(dFP[0]*dFP[0]+dFP[1]*dFP[1]))*dFP[1];
        //             // }
        //         }
        //         else
        //         {
        //             // here we take the gradient of the 2 neighbouring cell pressure.
        //             Cell *c1,*c2;
        //             c1=cn;
        //             c2=(e->neighbor[0]==cn)?(e->neighbor[1]):(e->neighbor[0]);
        //             float d12[2]={0.0f,0.0f};
        //             // d12 c1 --> c2
        //             for(int j=0;j<3;j++)
        //             {
        //                 d12[0]+=(c2->e[j]->vertex[0]->posn[0]+c2->e[j]->vertex[1]->posn[0])/6-(c1->e[j]->vertex[0]->posn[0]+c1->e[j]->vertex[1]->posn[0])/6;
        //                 d12[1]+=(c2->e[j]->vertex[0]->posn[1]+c2->e[j]->vertex[1]->posn[1])/6-(c1->e[j]->vertex[0]->posn[1]+c1->e[j]->vertex[1]->posn[1])/6;
        //             }

        //             e->vel[0]-=(dt/rho)*((c2->P-c1->P)/(d12[0]*d12[0]+d12[1]*d12[1]))*d12[0];
        //             e->vel[1]-=(dt/rho)*((c2->P-c1->P)/(d12[0]*d12[0]+d12[1]*d12[1]))*d12[1];

        //         }
        //     }
        // }

    }
    void find_performance()
    {
        float ForceV[2]={0.0,0.0};
        float ForceP[2]={0.0,0.0};
        for(int i=0;i<m.NUM_EDGE;i++)
        {
            Edge *e=&m.ME[i];
            if(e->isBounEdge && e->isInsideEdge)
            {
                Cell *cp;
                cp=(e->neighbor[0]==nullptr)?(e->neighbor[1]):(e->neighbor[0]);
                Eigen::Matrix3f mat_p;
                Eigen::Vector3f b_p_x;
                Eigen::Vector3f b_p_y;

                for(int i=0;i<3;i++)
                {
                    mat_p(i,0)=1;
                    mat_p(i,1)=(cp->e[i]->vertex[0]->posn[0]+cp->e[i]->vertex[1]->posn[0])/2;
                    mat_p(i,2)=(cp->e[i]->vertex[0]->posn[1]+cp->e[i]->vertex[1]->posn[1])/2;
                    b_p_x[i]=cp->e[i]->vel[0];
                    b_p_y[i]=cp->e[i]->vel[1];
                }

                Eigen::Vector3f p_x = mat_p.partialPivLu().solve(b_p_x);
                Eigen::Vector3f p_y = mat_p.partialPivLu().solve(b_p_y);

                float a_p[2],b_p[2],c_p[2];
                a_p[0]=p_x[0]; a_p[1]=p_y[0];
                b_p[0]=p_x[1]; b_p[1]=p_y[1];
                c_p[0]=p_x[2]; c_p[1]=p_y[2];

                float nf[2];
                Node *a=e->vertex[0];
                Node *b=e->vertex[1];
                Node *c;
                for(int i=0;i<3;i++)
                {
                    if(cp->vertex[i]!=a&&cp->vertex[i]!=b)
                    {
                        c=cp->vertex[i];
                        break;
                    }
                }
                find_normal(c->posn[0],c->posn[1],a->posn[0],a->posn[1],b->posn[0],b->posn[1],nf);
                nf[0]*=-1;
                nf[1]*=-1;
                float flen=sqrt(pow(a->posn[0]-b->posn[0],2)+pow(a->posn[1]-b->posn[1],2));
                ForceP[0]-=cp->P*nf[0]*flen;
                ForceP[1]-=cp->P*nf[1]*flen;
                float Tx=b_p[0]*nf[0]+c_p[0]*nf[1]+b_p[0]*nf[0]+b_p[1]*nf[1];
                float Ty=b_p[1]*nf[0]+c_p[1]*nf[1]+c_p[0]*nf[0]+c_p[1]*nf[1];
                ForceV[0]+=cp->visc*Tx*flen;
                ForceV[1]+=cp->visc*Ty*flen;
            }
        }
        cout<<ForceP[0]<<" "<<ForceP[1]<<endl;
        cout<<ForceV[0]<<" "<<ForceV[1]<<endl;
        Lift=ForceP[1]+ForceV[1];
        Drag=ForceP[0]+ForceV[0];
    }
};

void dumpEdges1(const Edge *edg, int NUM_EDGE, const char* filename)
{
    ofstream file(filename);

    if(!file.is_open())
    {
        cerr << "Could not open " << filename << endl;
        return;
    }

    int id = 0;

    for(int i = 0; i < NUM_EDGE; i++)
    {
        Node* A = edg[i].vertex[0];
        Node* B = edg[i].vertex[1];

        if(A == nullptr || B == nullptr)
        {
            cerr << "Warning: Edge " << i << " has null vertex!" << endl;
            continue;
        }

        file << id << " "
             << A->posn[0] << " "
             << A->posn[1] << " "
             << B->posn[0] << " "
             << B->posn[1] << " "
             << edg[i].isBounEdge << " "
             << edg[i].isPdef << " "
             << edg[i].isVdef
             << "\n";

        id++;
    }

    file.close();

    cout << "Edges dumped to " << filename << endl;
}

void dumpCFDResults(Cell *MC, int NUM_CELL,Edge *ME, int NUM_EDGE,string cname,string ename)
{
    // ============================================================
    // Cell pressure + triangle geometry
    // ============================================================

    ofstream cellFile(cname);

    if(!cellFile.is_open())
    {
        cerr << "Could not open cell_data.csv" << endl;
        return;
    }

    cellFile << "cell_id,"
             << "x1,y1,"
             << "x2,y2,"
             << "x3,y3,"
             << "xc,yc,"
             << "pressure\n";

    for(int i = 0; i < NUM_CELL; i++)
    {
        Node *n0 = MC[i].vertex[0];
        Node *n1 = MC[i].vertex[1];
        Node *n2 = MC[i].vertex[2];

        float x1 = n0->posn[0];
        float y1 = n0->posn[1];

        float x2 = n1->posn[0];
        float y2 = n1->posn[1];

        float x3 = n2->posn[0];
        float y3 = n2->posn[1];

        float xc = (x1 + x2 + x3) / 3.0f;
        float yc = (y1 + y2 + y3) / 3.0f;

        cellFile << i << ","
                 << x1 << "," << y1 << ","
                 << x2 << "," << y2 << ","
                 << x3 << "," << y3 << ","
                 << xc << "," << yc << ","
                 << MC[i].P << "\n";
    }

    cellFile.close();


    // ============================================================
    // Edge velocity
    // ============================================================

    ofstream edgeFile(ename);

    if(!edgeFile.is_open())
    {
        cerr << "Could not open edge_data.csv" << endl;
        return;
    }

    edgeFile << "edge_id,x,y,u,v,velocity_magnitude,is_boundary\n";

    for(int i = 0; i < NUM_EDGE; i++)
    {
        Node *A = ME[i].vertex[0];
        Node *B = ME[i].vertex[1];

        if(A == nullptr || B == nullptr)
            continue;

        float x = 0.5f * (A->posn[0] + B->posn[0]);
        float y = 0.5f * (A->posn[1] + B->posn[1]);

        float u = ME[i].vel[0];
        float v = ME[i].vel[1];

        float mag = sqrt(u*u + v*v);

        bool isBound=ME[i].isBounEdge;

        edgeFile << i << ","
                 << x << "," << y << ","
                 << u << "," << v << ","
                 << mag << ","
                 << isBound << "\n";
    }

    edgeFile.close();

    cout << "CFD results dumped successfully." << endl;
}

void savePoints(Node* nodes, int n)
{
    ofstream file("points.csv");

    if (!file.is_open())
    {
        cerr << "Could not open points.csv\n";
        return;
    }

    file << "x,y\n";

    for (int i = 0; i < n; i++)
    {
        file << nodes[i].posn[0] << ","
             << nodes[i].posn[1] << "\n";
    }

    file.close();

    cout << "Saved " << n << " points to points.csv\n";
}
