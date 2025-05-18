module menu;

import surface_properties;
import geometry_loader;
import spectral_solver;
import radiation_solver;
import conduction_solver;
import visualizer;

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <chrono>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "imgui.h"
#include "imgui_impl_glut.h"
#include "imgui_impl_opengl2.h"

enum class AppState { MENU, SIMULATION };
static AppState appState= AppState::MENU;

static const int g_maxSteps= 5000;
static int   g_stepCount=0;
static double g_alpha= 0.0005; 
static double g_kappa= 0.25;  
int g_dt= 300;                

extern RadiationModel g_model;
static int ui_method=1;

static void keyboardCB(unsigned char key, int x, int y){
    if(!ImGui::GetIO().WantCaptureKeyboard){
        if(key==27){
            std::cout<<"Esc pressed. Exiting.\n";
            exit(0);
        }
    }
    ImGui_ImplGLUT_KeyboardFunc(key,x,y);
}
static void keyboardUpCB(unsigned char key, int x, int y){
    ImGui_ImplGLUT_KeyboardUpFunc(key,x,y);
}
static void specialCB(int key, int x, int y){
    ImGui_ImplGLUT_SpecialFunc(key,x,y);
}
static void specialUpCB(int key, int x, int y){
    ImGui_ImplGLUT_SpecialUpFunc(key,x,y);
}
static void mouseCB(int button, int state, int x, int y){
    ImGui_ImplGLUT_MouseFunc(button,state,x,y);
    ImGuiIO &io= ImGui::GetIO();
    if(state==GLUT_DOWN && (button==3||button==4)){
        io.MouseWheel+= (button==3)? 1: -1;
        return;
    }
    if(io.WantCaptureMouse) return;

    if(button==GLUT_LEFT_BUTTON){
        if(state==GLUT_DOWN){
            g_arcballDragging=true;
            g_virtualX=x;
            g_virtualY=y;
        } else {
            g_arcballDragging=false;
        }
    }
    if(button==GLUT_RIGHT_BUTTON){
        g_rightDown= (state==GLUT_DOWN);
    }
    g_lastMouseX= x;
    g_lastMouseY= y;
}
static void motionCB(int x, int y){
    ImGui_ImplGLUT_MotionFunc(x,y);
    ImGuiIO &io= ImGui::GetIO();
    if(io.WantCaptureMouse) return;

    mouseMotionArcball(x,y,
                       g_arcballDragging,
                       g_rightDown,
                       g_lastMouseX,
                       g_lastMouseY,
                       g_virtualX,
                       g_virtualY,
                       g_zoom,
                       g_sceneQuat);
    glutPostRedisplay();
}
static void reshapeCB(int w, int h){
    windowWidth= w;
    windowHeight= h;
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60., double(w)/double(h), 0.1,100.);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    ImGuiIO &io= ImGui::GetIO();
    io.DisplaySize= ImVec2((float)w,(float)h);
}

static void displayCB(){
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGLUT_NewFrame();
    ImGui::NewFrame();

    if(appState==AppState::MENU){
        ImGui::SetNextWindowPos(ImVec2(0,0));
        ImGui::SetNextWindowSize(ImVec2((float)windowWidth, (float)windowHeight));
        ImGui::Begin("Simulation Settings");

        ImGui::BeginChild("MainChild", ImVec2(0,0), true);
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(12,8));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(12,12));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(20,20));

        const char* methods[]={"Monte Carlo","Discrete Ordinates","S2S"};
        ImGui::ListBox("Radiation Model", &ui_method, methods, IM_ARRAYSIZE(methods), 3);
        switch(ui_method){
            case 0: g_model= RadiationModel::MONTE_CARLO; break;
            case 1: g_model= RadiationModel::DISCRETE_ORDINATES; break;
            case 2: g_model= RadiationModel::S2S; break;
        }

        {
            double alphaStep=0.0001;
            double alphaStepFast=0.001;
            ImGui::InputScalar("Alpha (radiation step)", ImGuiDataType_Double,
                               &g_alpha,
                               &alphaStep, &alphaStepFast, "%.6f",
                               ImGuiInputTextFlags_CharsDecimal|ImGuiInputTextFlags_CharsScientific);
        }
        {
            double kappaStep=0.01;
            double kappaStepFast=0.1;
            ImGui::InputScalar("Kappa (conduction step)", ImGuiDataType_Double,
                               &g_kappa,
                               &kappaStep, &kappaStepFast, "%.6f",
                               ImGuiInputTextFlags_CharsDecimal|ImGuiInputTextFlags_CharsScientific);
        }
        ImGui::InputInt("Simulation Step (ms)", &g_dt, 10,10);

        if(ImGui::CollapsingHeader("Radiation Model Parameters", ImGuiTreeNodeFlags_DefaultOpen)){
            if(g_model==RadiationModel::MONTE_CARLO){
                ImGui::InputInt("Monte Carlo Samples", &monteCarloSamples, 100,100);
            } else if(g_model==RadiationModel::DISCRETE_ORDINATES){
                ImGui::InputInt("nTheta", &discreteOrdinates_nTheta,1,1);
                ImGui::InputInt("nPhi", &discreteOrdinates_nPhi,1,1);
            } else if(g_model==RadiationModel::S2S){
                ImGui::InputInt("S2S Subdivision", &s2sSubDiv,1,1);
            }
        }

        if(ImGui::CollapsingHeader("Band Edges", ImGuiTreeNodeFlags_DefaultOpen)){
            static std::vector<double> bandEdgesUI{0.1, 2.0, 5.0, 20.0};
            double step= 1.0;
            double fastStep= 1.0;
            for(size_t i=0; i< bandEdgesUI.size(); i++){
                ImGui::PushID((int)i);
                ImGui::InputScalar("##BandEdge", ImGuiDataType_Double,
                                   &bandEdgesUI[i], &step, &fastStep, "%.6f");
                ImGui::SameLine();
                ImGui::Text("10^-6 m");
                ImGui::SameLine();
                if(ImGui::Button("Remove")){
                    bandEdgesUI.erase(bandEdgesUI.begin()+ i);
                    ImGui::PopID();
                    break;
                }
                ImGui::PopID();
            }
            if(ImGui::Button("Add Band")){
                double newEdge= bandEdgesUI.empty()? 1.0: bandEdgesUI.back()+1.0;
                bandEdgesUI.push_back(newEdge);
            }
        }

        if(ImGui::CollapsingHeader("Add Object", ImGuiTreeNodeFlags_DefaultOpen)){
            static char stlFile[256]= "object_ascii.stl";
            static float pos[3]={0.f, 0.f, 0.f};
            static double initTemp=300.0;
            static int bcTypeIdx=0;
            static double bcValue=0.0;
            static double objectMass=1.0;
            static double objectCp=1000.0;
            static double reflectCoeff=0.5;
            static double absorbCoeff=0.5;

            ImGui::InputText("STL File##add", stlFile, IM_ARRAYSIZE(stlFile));
            float posStep=1.f;
            float posFastStep=1.f;
            ImGui::InputScalarN("Position", ImGuiDataType_Float, pos, 3,
                                &posStep, &posFastStep, "%.2f");

            double initTempStep=10.0, initTempFastStep=10.0;
            ImGui::InputScalar("Initial Temperature##add", ImGuiDataType_Double,
                               &initTemp,
                               &initTempStep, &initTempFastStep, "%.6f",
                               ImGuiInputTextFlags_CharsDecimal|ImGuiInputTextFlags_CharsScientific);

            const char* bcTypes[]={"1","2"};
            ImGui::Combo("Boundary Condition##add", &bcTypeIdx, bcTypes, IM_ARRAYSIZE(bcTypes));
            if(bcTypeIdx==1){
                double fluxStep=10.0;
                double fluxFastStep=10.0;
                ImGui::InputScalar("Heat Flux (W/m²)##add", ImGuiDataType_Double,
                                   &bcValue, &fluxStep, &fluxFastStep, "%.6f",
                                   ImGuiInputTextFlags_CharsDecimal|ImGuiInputTextFlags_CharsScientific);
            }

            double massStep=0.1, massFastStep=0.1;
            ImGui::InputScalar("Mass (kg)##add", ImGuiDataType_Double,
                               &objectMass,
                               &massStep, &massFastStep, "%.6f",
                               ImGuiInputTextFlags_CharsDecimal|ImGuiInputTextFlags_CharsScientific);
            double cpStep=10.0, cpFastStep=10.0;
            ImGui::InputScalar("Specific Heat (J/kg*K)##add", ImGuiDataType_Double,
                               &objectCp,
                               &cpStep, &cpFastStep, "%.6f",
                               ImGuiInputTextFlags_CharsDecimal|ImGuiInputTextFlags_CharsScientific);

            double reflStep=0.1, reflFastStep=0.1;
            ImGui::InputScalar("Reflection Coefficient##add", ImGuiDataType_Double,
                               &reflectCoeff,
                               &reflStep, &reflFastStep, "%.3f",
                               ImGuiInputTextFlags_CharsDecimal|ImGuiInputTextFlags_CharsScientific);

            double absorpStep=0.1, absorpFastStep=0.1;
            ImGui::InputScalar("Absorption Coefficient##add", ImGuiDataType_Double,
                               &absorbCoeff,
                               &absorpStep, &absorpFastStep, "%.3f",
                               ImGuiInputTextFlags_CharsDecimal|ImGuiInputTextFlags_CharsScientific);

            if(ImGui::Button("Add Object##button")){
                Object newObject;
                std::vector<Vector3> vertices;
                std::vector<Face> faces;
                if(!loadSTL(stlFile, vertices, faces)){
                    std::cerr<<"Failed to load STL file: "<<stlFile<<"\n";
                } else {
                    for(auto &v: vertices){
                        v.x+= pos[0];
                        v.y+= pos[1];
                        v.z+= pos[2];
                    }
                    newObject.vertices= std::move(vertices);
                    newObject.faces= std::move(faces);
                    for(auto &fc: newObject.faces){
                        fc.temp= initTemp;
                    }
                    double totalA=0.;
                    for(auto &fc: newObject.faces){
                        totalA+= fc.area;
                    }
                    newObject.totalArea= totalA;
                    newObject.mass= objectMass;
                    newObject.cp  = objectCp;
                    for(int b=0; b<3; b++){
                        newObject.reflectivity[b]= reflectCoeff;
                        newObject.emissivity[b]  = absorbCoeff;
                    }
                    newObject.boundaryCondition= (bcTypeIdx==0)?BoundaryConditionType::First:BoundaryConditionType::Second;
                    newObject.boundaryValue= bcValue;
                    newObject.stlFileName= stlFile;
                    newObject.position= {pos[0], pos[1], pos[2]};
                    newObject.initialTemperature= initTemp;

                    g_objects.push_back(std::move(newObject));
                    updateGlobalData();
                    buildFaceGeomCache();

                    strcpy(stlFile,"object_ascii.stl");
                    pos[0]=pos[1]=pos[2]=0.f;
                    initTemp=300.;
                    bcTypeIdx=0;
                    bcValue=0.;
                    objectMass=1.0;
                    objectCp=1000.0;
                    reflectCoeff=0.5;
                    absorbCoeff=0.5;
                }
            }
        }

        int objectToRemove=-1;
        if(ImGui::CollapsingHeader("Objects", ImGuiTreeNodeFlags_DefaultOpen)){
            for(int i=0; i<(int)g_objects.size(); i++){
                ImGui::PushID(i);
                if(ImGui::TreeNode(("Object "+std::to_string(i+1)).c_str())){
                    ImGui::BeginChild("ObjectInfo", ImVec2(0,0), true);
                    ImGui::BeginDisabled(true);

                    {
                        char buf[256];
                        snprintf(buf, sizeof(buf), "%s", g_objects[i].stlFileName.c_str());
                        ImGui::InputText("STL File", buf, sizeof(buf));
                    }
                    {
                        char posBuffer[128];
                        snprintf(posBuffer,sizeof(posBuffer),"%.2f, %.2f, %.2f",
                                 g_objects[i].position.x,
                                 g_objects[i].position.y,
                                 g_objects[i].position.z);
                        ImGui::InputText("Position", posBuffer, sizeof(posBuffer));
                    }
                    {
                        char tempBuffer[64];
                        snprintf(tempBuffer,sizeof(tempBuffer),"%.2f", g_objects[i].initialTemperature);
                        ImGui::InputText("Initial Temp", tempBuffer, sizeof(tempBuffer));
                    }
                    {
                        const char* bcItems[]={"1","2"};
                        ImGui::Text("Boundary Condition: %s", bcItems[(int)g_objects[i].boundaryCondition]);
                        if(g_objects[i].boundaryCondition==BoundaryConditionType::Second){
                            char bValueBuffer[64];
                            snprintf(bValueBuffer,sizeof(bValueBuffer),"%.6f", g_objects[i].boundaryValue);
                            ImGui::InputText("Heat Flux (W/m²)", bValueBuffer, sizeof(bValueBuffer));
                        }
                        char massBuf[64];
                        snprintf(massBuf,sizeof(massBuf),"%.6f", g_objects[i].mass);
                        ImGui::InputText("Mass (kg)", massBuf, sizeof(massBuf));
                        char cpBuf[64];
                        snprintf(cpBuf,sizeof(cpBuf),"%.6f", g_objects[i].cp);
                        ImGui::InputText("Specific Heat (J/kg*K)", cpBuf, sizeof(cpBuf));
                        char areaBuf[64];
                        snprintf(areaBuf,sizeof(areaBuf),"%.6f", g_objects[i].totalArea);
                        ImGui::InputText("Total Area (m^2)", areaBuf, sizeof(areaBuf));
                    }
                    {
                        char reflBuf[64];
                        snprintf(reflBuf,sizeof(reflBuf),"%.3f", g_objects[i].reflectivity[0]);
                        ImGui::InputText("Reflection Coefficient", reflBuf, sizeof(reflBuf));
                        char absorpBuf[64];
                        snprintf(absorpBuf,sizeof(absorpBuf),"%.3f", g_objects[i].emissivity[0]);
                        ImGui::InputText("Absorption Coefficient", absorpBuf, sizeof(absorpBuf));
                    }

                    ImGui::EndDisabled();
                    ImGui::EndChild();

                    if(ImGui::Button("Remove Object")){
                        objectToRemove=i;
                        ImGui::TreePop();
                        ImGui::PopID();
                        break;
                    }
                    ImGui::TreePop();
                }
                ImGui::PopID();
            }
            if(objectToRemove>=0){
                g_objects.erase(g_objects.begin()+ objectToRemove);
                updateGlobalData();
            }
        }

        if(ImGui::Button("Start Simulation", ImVec2(0,50))){
            if(g_objects.empty()){
                std::cerr<<"No objects defined!\n";
            } else {
                appState=AppState::SIMULATION;
                g_stepCount=0;
                simulationRunning=true;

                updateMinMaxTemps();
                switch(g_model) {
                    case RadiationModel::DISCRETE_ORDINATES:
                        computeDONodes(discreteOrdinates_nTheta, discreteOrdinates_nPhi);
                        break;
                    case RadiationModel::S2S:
                        computeViewFactors_S2S(s2sSubDiv);
                        computeTriInfos(s2sSubDiv);
                        break;
                }

                std::thread simThread([](){
                    while(simulationRunning && g_stepCount<g_maxSteps){
                        {
                            std::lock_guard<std::mutex> lock(simMutex);
                            radiationFull(g_alpha);
                            conductionOneStepAll(g_kappa, g_dt / 100);
                            applyBoundaryConditions(0.1);
                            g_stepCount++;
                        }
                        std::this_thread::sleep_for(std::chrono::milliseconds(g_dt));
                    }
                    simulationRunning=false;
                    std::cout<<"Simulation done, steps="<<g_stepCount<<"\n";
                });
                simThread.detach();
            }
        }

        ImGui::PopStyleVar(3);
        ImGui::EndChild();
        ImGui::End();
    }

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    if(appState==AppState::SIMULATION){
        std::lock_guard<std::mutex> lock(simMutex);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0,0,-10.f*g_zoom);

        double rotMatrix[16];
        {
            double xx=g_sceneQuat.x*g_sceneQuat.x, yy=g_sceneQuat.y*g_sceneQuat.y, zz=g_sceneQuat.z*g_sceneQuat.z;
            double xy=g_sceneQuat.x*g_sceneQuat.y, xz=g_sceneQuat.x*g_sceneQuat.z, yz=g_sceneQuat.y*g_sceneQuat.z;
            double wx=g_sceneQuat.w*g_sceneQuat.x, wy=g_sceneQuat.w*g_sceneQuat.y, wz=g_sceneQuat.w*g_sceneQuat.z;
            rotMatrix[0]= 1-2*(yy+zz); rotMatrix[1]= 2*(xy+wz);   rotMatrix[2]= 2*(xz-wy);   rotMatrix[3]=0;
            rotMatrix[4]= 2*(xy-wz);   rotMatrix[5]= 1-2*(xx+zz); rotMatrix[6]= 2*(yz+wx);   rotMatrix[7]=0;
            rotMatrix[8]= 2*(xz+wy);   rotMatrix[9]= 2*(yz-wx);   rotMatrix[10]=1-2*(xx+yy); rotMatrix[11]=0;
            rotMatrix[12]=0; rotMatrix[13]=0; rotMatrix[14]=0; rotMatrix[15]=1;
        }
        glMultMatrixd(rotMatrix);

        renderAllFaces(currentMinTemp, currentMaxTemp);
        renderColorBar(windowWidth, windowHeight, currentMinTemp, currentMaxTemp);
    }

    ImGui::Render();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
    glutSwapBuffers();
}

static void idleCB(){
    glutPostRedisplay();
}

int main(int argc, char** argv){
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_DEPTH|GLUT_RGBA);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("Radiation Simulation Example");

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io= ImGui::GetIO();
    io.Fonts->Clear();
    ImFont* customFont= io.Fonts->AddFontFromFileTTF("resources/Roboto-Medium.ttf",18.f);
    if(!customFont){
        std::cout<<"Failed to load custom font, using default.\n";
    }
    io.FontGlobalScale=1.0f;

    ImGui::StyleColorsLight();
    ImGuiStyle& style = ImGui::GetStyle();

    style.WindowRounding = 2.0f;
    style.FrameRounding  = 2.0f;
    style.GrabRounding   = 2.0f;

    style.Colors[ImGuiCol_WindowBg]   = ImVec4(0.95f, 0.95f, 0.95f, 1.0f);
    style.Colors[ImGuiCol_ChildBg]    = ImVec4(0.90f, 0.90f, 0.90f, 1.0f); 
    style.Colors[ImGuiCol_Border]     = ImVec4(0.80f, 0.80f, 0.80f, 1.0f); 
    style.Colors[ImGuiCol_Text]       = ImVec4(0.0f, 0.0f, 0.0f, 1.0f); 

    style.Colors[ImGuiCol_Button]         = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
    style.Colors[ImGuiCol_ButtonHovered]  = ImVec4(0.75f, 0.75f, 0.75f, 1.0f);
    style.Colors[ImGuiCol_ButtonActive]   = ImVec4(0.70f, 0.70f, 0.70f, 1.0f);
        
    ImGui_ImplGLUT_Init();
    ImGui_ImplOpenGL2_Init();

    initGL();

    glutKeyboardFunc(keyboardCB);
    glutKeyboardUpFunc(keyboardUpCB);
    glutSpecialFunc(specialCB);
    glutSpecialUpFunc(specialUpCB);
    glutDisplayFunc(displayCB);
    glutReshapeFunc(reshapeCB);
    glutMouseFunc(mouseCB);
    glutMotionFunc(motionCB);
    glutIdleFunc(idleCB);

    glutMainLoop();

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGLUT_Shutdown();
    ImGui::DestroyContext();
    if(g_bvhRoot) {
        deleteBVH(g_bvhRoot);
        g_bvhRoot = nullptr;
    }
    return 0;
}
