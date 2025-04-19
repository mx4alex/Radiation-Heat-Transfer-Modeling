APP = ThermalSimulation
CXX = g++
CXXFLAGS = -std=c++17 -fopenmp
IMGUI_DIR = C:/imgui
IMGUI_BACKENDS_DIR = C:/imgui/backends

PROJECT_SRCS = \
  geometry_loader.cpp \
  surface_properties.cpp \
  spectral_solver.cpp \
  radiation_solver.cpp \
  conduction_solver.cpp \
  visualizer.cpp \
  menu.cpp

IMGUI_SRCS = \
  $(IMGUI_DIR)/imgui.cpp \
  $(IMGUI_DIR)/imgui_draw.cpp \
  $(IMGUI_DIR)/imgui_widgets.cpp \
  $(IMGUI_DIR)/imgui_tables.cpp \
  $(IMGUI_BACKENDS_DIR)/imgui_impl_glut.cpp \
  $(IMGUI_BACKENDS_DIR)/imgui_impl_opengl2.cpp

SRCS = $(PROJECT_SRCS) $(IMGUI_SRCS)
LIBS = -lfreeglut -lopengl32 -lglu32
INCLUDES = -I$(IMGUI_DIR) -I$(IMGUI_BACKENDS_DIR)

$(APP): $(SRCS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRCS) -o $(APP) $(LIBS)

clean:
	rm -f $(APP)
