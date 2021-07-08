
#include "delfem2/glfw/viewer2.h"
#include "delfem2/mshprimitive.h"
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <cstdio>
#include <Eigen/Core>
#include <Eigen/SVD>
//#include <Eigen/LU>

/**
 * making 2D regular grid with edge length 1
 * @tparam REAL float or double
 * @param aXYZ list of points (XYXYXY...)
 * @param aQuad vertex index for quad mesh
 * @param nx nubmer of quads in x direction
 * @param ny number of quads in y direction
 */
template <typename REAL>
void MeshQuad2_Grid(
    std::vector<REAL>& aXYZ,
    std::vector<unsigned int>& aQuad,
    unsigned int nx,
    unsigned int ny)
{
  unsigned int np = (nx+1)*(ny+1);
  aXYZ.resize(np*2);
  for(unsigned int iy=0;iy<ny+1;++iy){
    for(unsigned int ix=0;ix<nx+1;++ix){
      unsigned int ip = iy*(nx+1)+ix;
      aXYZ[ip*2+0] = ix;
      aXYZ[ip*2+1] = iy;
    }
  }
  aQuad.resize(nx*ny*4);
  for(unsigned int iy=0;iy<ny;++iy){
    for(unsigned int ix=0;ix<nx;++ix){
      unsigned int iq = iy*nx+ix;
      aQuad[iq*4+0] = (iy+0)*(nx+1)+(ix+0);
      aQuad[iq*4+1] = (iy+0)*(nx+1)+(ix+1);
      aQuad[iq*4+2] = (iy+1)*(nx+1)+(ix+1);
      aQuad[iq*4+3] = (iy+1)*(nx+1)+(ix+0);
    }
  }
}

// ----------------------------------

// print out error
static void error_callback( int error, const char* description){ fputs(description, stderr); }

int main()
{
  std::vector<float> aXY; // coordinates
  std::vector<unsigned int> aQuad; // index of points
  std::vector<float> aMass; // point mass
  { // initialize mesh and point masses
    const unsigned int nx = 3;
    const unsigned int ny = 10;
    const float scale = 0.1;
    // make grid
    MeshQuad2_Grid(
        aXY, aQuad,
        3, 10);
    // set mass
    aMass.resize(aXY.size() / 2, 1.f);
    for (unsigned int ixy = 0; ixy < aXY.size() / 2; ++ixy){
      if( aXY[ixy*2+1] > 0.001 ){ continue; }
      aMass[ixy] = 10000.0; // mass at the boundary is very large to prevent moving during shape matching operation
    }
    for(unsigned int ixy=0;ixy<aXY.size()/2;++ixy) {
      aXY[ixy*2+0] = (aXY[ixy*2+0] - nx*0.5)*scale;
      aXY[ixy*2+1] = (aXY[ixy*2+1] - ny*0.5)*scale;
    }
  }
  const std::vector<float> aXY0 = aXY; // rest shape
  std::vector<float> aUV(aXY.size(),0.0); // velocity
  std::vector<float> aXYt(aXY.size(),0.0); // tentative shape
  std::vector<float> aXYt_update = aXY;
  std::vector<float> aNO(aQuad.size(), 0.0);
  float dt = 0.03;
  // ---------
  delfem2::glfw::CViewer2 viewer;
  {
    viewer.view_height = 1.0;
    viewer.title = "task9: Shape Matching Deformation";
  }
  glfwSetErrorCallback(error_callback);
  if ( !glfwInit() ) { exit(EXIT_FAILURE); }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  // ------
  viewer.InitGL();
  unsigned int c = 0;
  while (!glfwWindowShouldClose(viewer.window))
  {
    { // animation using Shape Matching Deformation [Müller et al. 2005]
      for(unsigned int ixy=0;ixy<aXY.size()/2;++ixy) { // set position & velocity at the boundary
        if( aMass[ixy] < 1.1 ){ continue; }
        aXY[ixy*2+0] = aXY0[ixy*2+0]+0.2*sin(glfwGetTime()*5);
        aXY[ixy*2+1] = aXY0[ixy*2+1];
        aUV[ixy*2+0] = 0.0;
        aUV[ixy*2+1] = 0.0;
      }
      for(unsigned int ixy=0;ixy<aXY.size()/2;++ixy){ // compute tentative shape
        aXYt[ixy*2+0] = aXY[ixy*2+0] + dt*aUV[ixy*2+0];
        aXYt[ixy*2+1] = aXY[ixy*2+1] + dt*aUV[ixy*2+1];
      }
        std::vector<int> flag(aXY.size(),0);
        std::cout << "aXY.size() = " << aXY.size() << std::endl;
        std::cout << "aQuad.size() = " << aQuad.size() << std::endl;

        for(unsigned k=0;k<aXY.size();++k){
          std::cout << " f=" << flag[k];
        }
        std::cout << std::endl;

        // // print out aMass
        // for(unsigned k=0;k<aXY.size()/2;++k){
        //   std::cout << "am" << k << "=" << aMass[k] << std::endl;
        // }

        // // print out aUV
        // for(unsigned int k=0;k<aUV.size();++k){
        //   std::cout << "aUV" << k << "=" << aUV[k] << std::endl;
        // }

      for(unsigned int iq=0;iq<aQuad.size()/4;++iq){
        const Eigen::Vector2f ap[4] = { // coordinates of quad's corner points (tentative shape)
            Eigen::Map<Eigen::Vector2f>(aXYt.data()+aQuad[iq*4+0]*2),
            Eigen::Map<Eigen::Vector2f>(aXYt.data()+aQuad[iq*4+1]*2),
            Eigen::Map<Eigen::Vector2f>(aXYt.data()+aQuad[iq*4+2]*2),
            Eigen::Map<Eigen::Vector2f>(aXYt.data()+aQuad[iq*4+3]*2) };
        const Eigen::Vector2f aq[4] = { // coordinates of quads' corner points (rest shape)
            Eigen::Map<const Eigen::Vector2f>(aXY0.data()+aQuad[iq*4+0]*2),
            Eigen::Map<const Eigen::Vector2f>(aXY0.data()+aQuad[iq*4+1]*2),
            Eigen::Map<const Eigen::Vector2f>(aXY0.data()+aQuad[iq*4+2]*2),
            Eigen::Map<const Eigen::Vector2f>(aXY0.data()+aQuad[iq*4+3]*2) };
        const float am[4] = { // masses of the points
            aMass[aQuad[iq*4+0]],
            aMass[aQuad[iq*4+1]],
            aMass[aQuad[iq*4+2]],
            aMass[aQuad[iq*4+3]] };
        // write some code below to rigidly transform the points in the rest shape (`aq`) such that the
        // weighted sum of squared distances against the points in the tentative shape (`ap`) is minimized (`am` is the weight).
        Eigen::Vector2f tcg = (am[0]*ap[0]+am[1]*ap[1]+am[2]*ap[2]+am[3]*ap[3])/(am[0]+am[1]+am[2]+am[3]);
        Eigen::Vector2f Tcg = (am[0]*aq[0]+am[1]*aq[1]+am[2]*aq[2]+am[3]*aq[3])/(am[0]+am[1]+am[2]+am[3]);
        auto BAT = am[0]*(ap[0]-tcg)*(aq[0]-Tcg).transpose()+am[1]*(ap[1]-tcg)*(aq[1]-Tcg).transpose()+am[2]*(ap[2]-tcg)*(aq[2]-Tcg).transpose()+am[3]*(ap[3]-tcg)*(aq[3]-Tcg).transpose();
        
        // svd
        Eigen::JacobiSVD<Eigen::MatrixXf> svd(BAT , Eigen::ComputeThinV | Eigen::ComputeThinU );
        auto U = svd.matrixU();
        auto V = svd.matrixV();

        //Eigen::MatrixXf Ropt = U*V.transpose();
        Eigen::MatrixXf Ropt = U*V.transpose();
        Eigen::Vector2f topt = tcg - Ropt*Tcg;
        Eigen::Vector2f update(topt.size(),0.0);

        // for(unsigned int j=0;j<4;++j){
        //   update = Ropt*aq[j] - topt;
        //   aXYt[aQuad[iq*4+j]*2] = update[0];
        //   aXYt[aQuad[iq*4+j]*2+1] = update[1];
        // }


        //update 
        // for(unsigned int j=0;j<4;++j){
        //   //if(flag[aQuad[iq*4+j]*2] == 1){continue;}
        //   Eigen::Vector2f Xupdate = Ropt*aq[j] - topt;
        //   if(aXYt[aQuad[iq*4+j]*2] != Xupdate[0]){aXYt[aQuad[iq*4+j]*2] = Xupdate[0];flag[aQuad[iq*4+j]*2] += 1;}
        //   if(aXYt[aQuad[iq*4+j]*2+1] != Xupdate[1]){aXYt[aQuad[iq*4+j]*2+1] = Xupdate[1];flag[aQuad[iq*4+j]*2+1] += 1;}
        //   // flag[aQuad[iq*4+j]*2] += 1;
        //   // flag[aQuad[iq*4+j]*2+1] += 1;
        // }

        
        // for(unsigned int j=0;j<4;j++){
        //   Eigen::Vector2f Xupdate = Ropt*aq[j] - topt;
        //   if(flag[aQuad[iq*4+j]*2]>=1){continue;}
        //   aXYt_update[aQuad[iq*4+j]*2] = Xupdate[0];
        //   aXYt_update[aQuad[iq*4+j]*2+1] = Xupdate[1];
        //   flag[aQuad[iq*4+j]*2] += 1;
        //   flag[aQuad[iq*4+j]*2+1] += 1;          
        // }
        
      //  for(unsigned int j=0;j<4;j++){
      //     Eigen::Vector2f Xupdate = Ropt*aq[j] - topt;
      //     //if(flag[aQuad[iq*4+j]*2]>=1){continue;}
      //     aXYt[aQuad[iq*4+j]*2] = Xupdate[0];
      //     aXYt[aQuad[iq*4+j]*2+1] = Xupdate[1];
      //     flag[aQuad[iq*4+j]*2] += 1;
      //     flag[aQuad[iq*4+j]*2+1] += 1;          
      //   }


        std::cout << iq << std::endl;
        for (unsigned int j=0;j<flag.size();j++)
        {
          std::cout << " f=" << flag[j];
        }
        std::cout << std::endl;

        // std:: cout << " = " << std::endl << ap[0] << std::endl;
        // std:: cout << "iq = " << std::endl << iq << std::endl;
        // std:: cout << "aQuad[iq*4+0]*2 " << std::endl << aQuad[iq*4+0]*2 << std::endl;
        // std:: cout << "(aXY0.data()+aQuad[iq*4+0]*2) " << std::endl << Eigen::Map<Eigen::Vector2f>(aXYt.data()+aQuad[iq*4+0]*2) << std::endl;
        // std:: cout << "ap = " << std::endl << ap[0] << std::endl;
        // std::cout << "axyt x" << std::endl << aXYt[aQuad[iq*4+0]*2] << std::endl;
        // std::cout << "axyt y" << std::endl << aXYt[aQuad[iq*4+0]*2+1] << std::endl;



        
        // no edits further down
      }
      //aXYt = aXYt_update;
      // aXY0 = aXYt;
      for(unsigned int ixy=0;ixy<aXY.size()/2;++ixy) { // update position and velocities
        aUV[ixy*2+0] = (aXYt[ixy*2+0] - aXY[ixy*2+0])/dt;
        aUV[ixy*2+1] = (aXYt[ixy*2+1] - aXY[ixy*2+1])/dt;
        aXY[ixy*2+0] = aXYt[ixy*2+0];
        aXY[ixy*2+1] = aXYt[ixy*2+1];
      }
    }
    // ----------------
    viewer.DrawBegin_oldGL();
    ::glColor3d(0,0,0);
    ::glBegin(GL_LINES);
    for(unsigned int iq=0;iq<aQuad.size()/4;++iq){ // draw edges of quads
      if(iq==0){::glColor3d(1,0,0);}
      ::glVertex2fv(aXY.data()+aQuad[iq*4+0]*2);
      ::glVertex2fv(aXY.data()+aQuad[iq*4+1]*2);
      ::glVertex2fv(aXY.data()+aQuad[iq*4+1]*2);
      ::glVertex2fv(aXY.data()+aQuad[iq*4+2]*2);
      ::glVertex2fv(aXY.data()+aQuad[iq*4+2]*2);
      ::glVertex2fv(aXY.data()+aQuad[iq*4+3]*2);
      ::glVertex2fv(aXY.data()+aQuad[iq*4+3]*2);
      ::glVertex2fv(aXY.data()+aQuad[iq*4+0]*2);
      ::glColor3d(0,0,0);

     

      // ::glVertex2fv(aNO.data()+iq*4+0);
      // ::glVertex2fv(aNO.data()+iq*4+1);
      // ::glVertex2fv(aNO.data()+iq*4+1);
      // ::glVertex2fv(aNO.data()+iq*4+2);
      // ::glVertex2fv(aNO.data()+iq*4+2);
      // ::glVertex2fv(aNO.data()+iq*4+3);
      // ::glVertex2fv(aNO.data()+iq*4+3);
      // ::glVertex2fv(aNO.data()+iq*4+0);
    }
    ::glEnd();
    viewer.SwapBuffers();
    glfwPollEvents();

    // for(unsigned int k=0;k<aXY.size()/2;++k){
    //   std::cout << "x = " << aXY[k*2+0] << "  y = " << aXY[k*2+1] << std::endl;
    // }
    // for(unsigned int k=0;k<aXY.size()/2;++k){
    //   std::cout << "xt = " << aXYt[k*2+0] << "  yt = " << aXYt[k*2+1] << std::endl;
    // }
    // for(unsigned int k=0;k<aXY.size()/2;++k){
    //   std::cout << "vx = " << aUV[k*2+0] << "  vy = " << aUV[k*2+1] << std::endl;
    // }

    // c++;
    // if(c==20){while(1);}
    
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}
