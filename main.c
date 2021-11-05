#include "global.h"
#include <cairo.h>
#include <gtk/gtk.h>
#include <math.h>
#include "messages.h"
#include "genes.h"
#include "array.h"

struct RCPair{
  int row;
  int col;
};

typedef Array(struct RCPair*) Array_RCPair;

Array_RCPair rcArray;

GtkBuilder *builder;

GtkWidget *window;

GtkWidget *grid_probabilities;
GtkWidget *lbl_probabilitiesPages;
int grid_probabilitiesRows = 0;
int lastRow = 0, lastCol = 0;
int desActualPageX = 1, desActualPageY = 1;
int desTotPageX = 1, desTotPageY = 1;
int maxRows = 5, maxCols = 5;
int totalProbs = 0;
bool needFillGrid = true;
double x = 0.0, y = 0.0;
double beginX = 0.0, beginY = 0.0;
bool pressed = false, geneNameAutoChanged = false;;

GtkWidget *btn_nFeatures;
GtkWidget *txt_nFeatures;
GtkWidget *txt_fileName;

GtkWidget *btn_loadFile;
GtkWidget *btn_saveFile;
GtkWidget *btn_chooseFile;
bool fileChosen = false;

GtkWidget *da_chromosomes;
float scale = 1.0;

bool autoFilledTxtGrid = false, autoFilledTxtFeatures = false;

void do_drawing(cairo_t *cr);
void createProbabilitiesGrid();
void fillProbabilitiesGrid(int total);
void fillProbabilitiesGridAux();
void initRCPairs();
void changeCommaToPoint(char str[9]);
void setEntriesRestrictions();

void on_txt_probability_changed(GtkEditable* self, gpointer user_data);
void on_txt_probability_insert_text(GtkEditable* self, gchar* new_text, gint new_text_length, gint* position, gpointer user_data);
void on_txt_gene_name_changed(GtkEditable* self, gpointer user_data);
void on_txt_probability_delete_text(GtkEditable* self,gint start_pos,gint end_pos, gpointer user_data);

int main(int argc, char *argv[]){
  initRCPairs();
  gtk_init(&argc, &argv);

  builder = gtk_builder_new();
  gtk_builder_add_from_file(builder, "glade/window.glade", NULL);

  window = GTK_WIDGET(gtk_builder_get_object(builder, "window"));
  grid_probabilities = GTK_WIDGET(gtk_builder_get_object(builder, "grid_probabilities"));
  da_chromosomes = GTK_WIDGET(gtk_builder_get_object(builder, "da_chromosomes"));
  txt_nFeatures = GTK_WIDGET(gtk_builder_get_object(builder, "txt_nFeatures"));
  txt_fileName = GTK_WIDGET(gtk_builder_get_object(builder, "txt_fileName"));
  lbl_probabilitiesPages = GTK_WIDGET(gtk_builder_get_object(builder, "lbl_probabilitiesPages"));

  gtk_window_maximize(GTK_WINDOW(window));
  gtk_builder_connect_signals(builder, NULL);
  createProbabilitiesGrid();

  setEntriesRestrictions();
  
  gtk_widget_show(window);
  gtk_main();
  return 0;
}

void setEntriesRestrictions(){
  gtk_entry_set_max_length (GTK_ENTRY(txt_nFeatures), 1);
}

void initRCPairs(){
  initArray(rcArray, struct RCPair*, maxRows*maxCols);
  for(int i = 0; i < rcArray.size; i++){
    rcArray.data[i] = NULL;
  }
}

gboolean on_da_chromosomes_draw(GtkWidget *widget, cairo_t *cr, gpointer user_data){
  do_drawing(cr);

  return FALSE;
}

void do_drawing(cairo_t *cr){
  float xCoordChromosome, yCoordGene;
  Array_char geneName;
  if(areGenesInitialized() && areOrdersCreated()){
    Array_ints orders = getOrders();
    cairo_scale(cr, scale, scale);
    cairo_set_line_width(cr, 1);
    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
    for(int i = 0; i < getTotalPossibleOrders(); i++){
      yCoordGene = 20+y;
      xCoordChromosome = 25+i*80+x;
      Array_float distances = getDistances(i);
      cairo_move_to(cr, xCoordChromosome, 10+y);
      cairo_line_to(cr, xCoordChromosome, 310+y);

      cairo_move_to(cr, xCoordChromosome-10, yCoordGene);
      cairo_line_to(cr, xCoordChromosome+10, yCoordGene);
      cairo_move_to(cr, xCoordChromosome+15, yCoordGene);
      geneName = getGeneName(orders.data[i].data[0]);
      cairo_show_text(cr, geneName.data);
      for(int j = 0; j < distances.used; j++){
        yCoordGene = yCoordGene+600*distances.data[j];
        cairo_move_to(cr, xCoordChromosome-10, yCoordGene);
        cairo_line_to(cr, xCoordChromosome+10, yCoordGene);
        cairo_move_to(cr, xCoordChromosome+15, yCoordGene);
        geneName = getGeneName(orders.data[i].data[j+1]);
        cairo_show_text(cr, geneName.data);
      }
    }
    cairo_stroke(cr);
  }
}

void createProbabilitiesGrid(){
  for(int i = 0; i <= maxRows; i++){
    gtk_grid_insert_row(GTK_GRID(grid_probabilities), i);
  }
  for(int i = 0; i <= maxCols; i++){
    gtk_grid_insert_column(GTK_GRID(grid_probabilities), i);
  }
  for(int row = 0; row <= maxRows; row++){
    for(int col = 0; col <= maxCols; col++){
      GtkWidget* cell = gtk_grid_get_child_at(GTK_GRID(grid_probabilities), col, row);
      if(!cell){
        cell = gtk_entry_new();
        gtk_grid_attach(GTK_GRID(grid_probabilities), cell, col, row, 1, 1);
        gtk_widget_set_sensitive(cell, FALSE);
        gtk_widget_show(cell);
        struct RCPair *rcp = malloc(sizeof(struct RCPair));
        rcp->row = row-1;
        rcp->col = col-1;
        appendPArray(rcArray, struct RCPair*, rcp);
        if(row != 0 && col != 0){
          g_signal_connect(cell, "insert-text", G_CALLBACK(on_txt_probability_insert_text), rcArray.data[rcArray.used-1]);
          g_signal_connect(cell, "changed", G_CALLBACK(on_txt_probability_changed), rcArray.data[rcArray.used-1]);
        } else {
          g_signal_connect(cell, "changed", G_CALLBACK(on_txt_gene_name_changed), rcArray.data[rcArray.used-1]);
          gtk_entry_set_max_length (GTK_ENTRY(cell), 7);
        }
      }
    }
  }
  gtk_grid_set_row_spacing(GTK_GRID(grid_probabilities), 10);
  gtk_grid_set_column_spacing(GTK_GRID(grid_probabilities), 10);
}

void fillProbabilitiesGrid(int total){
  if(totalProbs != total || needFillGrid){
    x = 0.0, y = 0.0;
    needFillGrid = false;
    grid_probabilitiesRows = totalProbs = total;
    lastRow = 0; lastCol = 0;
    desActualPageX = 1; desActualPageY = 1;
    desTotPageX = 1; desTotPageY = 1;
    fillProbabilitiesGridAux();
  }
}

void fillProbabilitiesGridAux(){
  int lRow = lastRow;
  int lCol = lastCol;
  int probRow=0,probCol=0;
  float prob = 0.0;
  char fillWith[9];
  geneNameAutoChanged = true;
  for(int row = 0; row <= maxRows; row++){
    lRow = row + lastRow;
    probRow = row+lastRow-1;
    for(int col = 0; col <= maxCols; col++){
      GtkWidget* cell = gtk_grid_get_child_at(GTK_GRID(grid_probabilities), col, row);
      lCol = col + lastCol;
      probCol = col+lastCol-1;
      if(row == 0 && col > 0){
        if(lCol-1 < getTotalGenes()){
          Array_char geneName = getGeneName(lCol-1);
          gtk_entry_set_text(GTK_ENTRY(cell), geneName.data);
          gtk_widget_set_sensitive(cell, TRUE);
        } else {
          gtk_entry_set_text(GTK_ENTRY(cell), "");
          gtk_widget_set_sensitive(cell, FALSE);
        }
      } else if(col == 0 && row > 0){
        if(lRow-1 < getTotalGenes()){
          Array_char geneName = getGeneName(lRow-1);
          gtk_entry_set_text(GTK_ENTRY(cell), geneName.data);
          gtk_widget_set_sensitive(cell, TRUE);
        } else {
          gtk_entry_set_text(GTK_ENTRY(cell), "");
          gtk_widget_set_sensitive(cell, FALSE);
        }
      } else if(row > 0 && col > 0){
        if(lCol <= lRow){
          gtk_entry_set_text(GTK_ENTRY(cell), "");
          gtk_widget_set_sensitive(cell, FALSE);
        } else {
          if(probRow < getTotalGenes() && probCol < getTotalGenes()){
            prob = getProbability(probRow, probCol);
            if(prob > 0 && prob < 0.5 && prob != 0.0 && prob != -0.0){
              autoFilledTxtGrid = true;
              sprintf(fillWith, "%f", prob);
              bye0sRight(fillWith);
              changeCommaToPoint(fillWith);
              gtk_entry_set_text(GTK_ENTRY(cell), fillWith);
              autoFilledTxtGrid = false;
            } else {
              gtk_entry_set_text(GTK_ENTRY(cell), "");
            }
            gtk_widget_set_sensitive(cell, TRUE);
          } else {
            gtk_entry_set_text(GTK_ENTRY(cell), "");
            gtk_widget_set_sensitive(cell, FALSE);
          }
        }
      } else {
        gtk_widget_set_sensitive(cell, FALSE);
      }
    }
  }
  geneNameAutoChanged = false;
  lastRow = lRow;
  lastCol = lCol;
  char lbl[23];
  if(grid_probabilitiesRows % maxRows == 0){
    desTotPageY = grid_probabilitiesRows / maxRows;
  } else{
    desTotPageX = grid_probabilitiesRows / maxRows + 1;
  }
  if(grid_probabilitiesRows % maxCols == 0){
    desTotPageY = grid_probabilitiesRows / maxCols;
  } else{
    desTotPageY = grid_probabilitiesRows / maxCols + 1;
  }
  sprintf(lbl, "(%d,%d) de (%d,%d)", desActualPageX, desActualPageY, desTotPageX, desTotPageY);
  gtk_label_set_text(GTK_LABEL(lbl_probabilitiesPages), lbl);
}

void on_window_destroy(){
  if(areGenesInitialized()){
    freeGenes();
  }
  for(int i = 0; i < rcArray.size; i++){
    if(rcArray.data[i]){
      free(rcArray.data[i]);
    }
  }
  freeArray(rcArray);
  g_object_unref(builder);
  gtk_main_quit();
}

void on_btn_chooseFile_file_set(GtkFileChooserButton *f){
  gchar *fileName = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(f));
  if(!fileChosen){
    fileChosen = true;
    if(!loadFile(fileName)){
      messagesWindow(window, "Could not open file");
    } else {
      needFillGrid = true;
      x = 0.0;
      y = 0.0;
      char fillWith[9];
      sprintf(fillWith, "%d", getTotalGenes());
      gtk_entry_set_text(GTK_ENTRY(txt_nFeatures), fillWith);
      fillProbabilitiesGrid(getTotalGenes());
      gtk_widget_queue_draw(window);
    }
    fileChosen = false;
  }
}

void on_txt_nFeatures_insert_text(GtkEditable* self, gchar* new_text, gint new_text_length, int position, gpointer user_data){
  if(!autoFilledTxtFeatures){
    if(new_text[new_text_length-1] < '2' || new_text[new_text_length-1] > '7'){
      g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
    }
  }
}

void on_btn_nFeatures_clicked(GtkButton *b){
  const gchar* val = gtk_entry_get_text(GTK_ENTRY(txt_nFeatures));
  int n = 0;
  if(gtk_entry_get_text_length(GTK_ENTRY(txt_nFeatures)) > 0){
    sscanf(val, "%d", &n);
  }
  if(n > 0){
    initializeGenes(n);
    fillProbabilitiesGrid(n);
  } else {
    messagesWindow(window, "Number of genes can not be empty");
  }
}

void on_btn_up_clicked(GtkButton *b){
  if(lastRow - 2*maxRows >= 0){
    lastRow = lastRow - 2*maxRows;
    lastCol = lastCol - maxCols;
    desActualPageX--;
    fillProbabilitiesGridAux();
  }
}

void on_btn_left_clicked(GtkButton *b){
  if(lastCol - 2*maxCols >= 0){
    lastRow = lastRow - maxRows;
    lastCol = lastCol - 2*maxCols;
    desActualPageY--;
    fillProbabilitiesGridAux();
  }
}

void on_btn_down_clicked(GtkButton *b){
  if(getTotalGenes() - lastRow > 0){
    lastCol = lastCol - maxCols;
    desActualPageX++;
    fillProbabilitiesGridAux();
  }
}

void on_btn_right_clicked(GtkButton *b){
  if(getTotalGenes() - lastCol > 0){
    lastRow = lastRow - maxRows;
    desActualPageY++;
    fillProbabilitiesGridAux();
  }
}

void needsChangePointToComma(char str[], int size){
  float a;
  sscanf("0.3", "%f", &a);
  if(a != 0.3){
    for(int i = 0; i < size; i++){
      if(str[i] == '.'){
        str[i] = ',';
      }
    }
  }
}

void changeCommaToPoint(char str[9]){
  for(int i = 0; i < 9; i++){
    if(str[i] == ','){
      str[i] = '.';
    }
  }
}

void on_txt_probability_changed(GtkEditable* self, gpointer user_data){
  float prob = 0.0;
  int row, col;
  char result[9];
  if(areGenesInitialized() && !autoFilledTxtGrid){
    const gchar* c = gtk_entry_get_text(GTK_ENTRY(self));
    strcpy(result,c);
    needsChangePointToComma(result, 9);
    sscanf(result, "%f", &prob);
    row = (*(struct RCPair*)user_data).row;
    col = (*(struct RCPair*)user_data).col;
    if(prob > 0.0 && prob < 0.5){
      insertProbability(prob, row + lastRow - maxRows, col + lastCol - maxCols);
    } else {
      insertProbability(2.0, row + lastRow - maxRows, col + lastCol - maxCols);
    }
    gtk_widget_queue_draw(window);
  }
}

void on_txt_probability_insert_text(GtkEditable* self, gchar* new_text, gint new_text_length, gint* position, gpointer user_data){
  float prob = 0.0;
  char result[9];
  if(!autoFilledTxtGrid){
    if(new_text_length <= 1){
      if(*position == 0){
        if(new_text[0] != '0'){
          g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
        }
      } else if(*position == 1){
        if(new_text[0] != '.'){
          g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
        }
      } else if(*position == 2){
        if(new_text[0] < '0' || new_text[0] >= '5'){
          g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
        }
      } else if(new_text[0] < '0' || new_text[0] > '9'){
        g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
      }
    } else {
      strcpy(result,new_text);
      needsChangePointToComma(result, 9);
      sscanf(result, "%f", &prob);
      if(prob <= 0.0 || prob >= 0.5){
        g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
      }
    }
  }
}

void  on_txt_gene_name_changed(GtkEditable* self, gpointer user_data){
  if(!geneNameAutoChanged){
    int row, col;
    row = (*(struct RCPair*)user_data).row;
    col = (*(struct RCPair*)user_data).col;
    printf("ROW:%d, COL:%d\n",row+1, col+1);
    const gchar* c = gtk_entry_get_text(GTK_ENTRY(self));
    char name[10];
    strcpy(name,c);
    if(row < 0){
      setGeneName(col, name);
    } else {
      setGeneName(row, name);
    }
    geneNameAutoChanged = true;
    GtkWidget* cell = gtk_grid_get_child_at(GTK_GRID(grid_probabilities), row+1, col+1);
    gtk_entry_set_text(GTK_ENTRY(cell), name);
    geneNameAutoChanged = false;
    gtk_widget_queue_draw(window);
  }
}

void on_btn_zoomIn_clicked(GtkButton *b){
  scale += 0.1;
  gtk_widget_queue_draw(window);
}

void on_btn_zoomOut_clicked(GtkButton *b){
  scale -= 0.1;
  gtk_widget_queue_draw(window);
}

void on_btn_inferProbabilities_clicked(GtkButton *b){
  if(inferProbabilities()){
    needFillGrid = true;
    fillProbabilitiesGrid(getTotalGenes());
    gtk_widget_queue_draw(window);
  } else {
    messagesWindow(window, "Could not infer probabilities");
  }
}

void on_btn_saveFile_clicked(GtkButton *b){
  const gchar* fn = gtk_entry_get_text(GTK_ENTRY(txt_fileName));
  bool isAlphaN = isalnum(fn[0]);
  if(isAlphaN && areGenesInitialized()){
    gchar* fileName = g_strconcat(fn, ".txt", NULL);
    if(!saveFile(fileName)){
      messagesWindow(window, "Could not save file");
    } else {
      messagesWindow(window, "File saved");
    }
  }
}

gboolean on_da_chromosomes_button_press_event(GtkWidget* self, GdkEventButton* event, gpointer user_data){
  beginX = event->x;
  beginY = event->y;
  pressed = true;
  return TRUE;
}

gboolean on_da_chromosomes_button_release_event(GtkWidget* self, GdkEventButton* event, gpointer user_data){
  pressed = false;
  return TRUE;
}

gboolean on_da_chromosomes_motion_notify_event(GtkWidget* self, GdkEventMotion *event, gpointer user_data){
  if(pressed){

  x += (event->x - beginX)/40;
  y += (event->y - beginY)/40;
  gtk_widget_queue_draw(window);
  }
  return TRUE;
}