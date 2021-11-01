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

GtkWidget *btn_nFeatures;
GtkWidget *txt_nFeatures;

GtkWidget *btn_loadFile;
GtkWidget *btn_saveFile;
GtkWidget *btn_chooseFile;
bool fileChosen = false;

GtkWidget *da_chromosomes;
double scale = 1.0;

bool needsDrawing = false;

void do_drawing(cairo_t *cr);
void createProbabilitiesGrid();
void fillProbabilitiesGrid(int total);
void fillProbabilitiesGridAux();
void initRCPairs();

void on_txt_probability_insert_text(GtkEditable* self, gchar* new_text, gint new_text_length, gint* position, gpointer user_data);


int main(int argc, char *argv[]){
  initRCPairs();
  gtk_init(&argc, &argv);

  builder = gtk_builder_new();
  gtk_builder_add_from_file(builder, "glade/window.glade", NULL);

  window = GTK_WIDGET(gtk_builder_get_object(builder, "window"));
  grid_probabilities = GTK_WIDGET(gtk_builder_get_object(builder, "grid_probabilities"));
  da_chromosomes = GTK_WIDGET(gtk_builder_get_object(builder, "da_chromosomes"));
  txt_nFeatures = GTK_WIDGET(gtk_builder_get_object(builder, "txt_nFeatures"));
  lbl_probabilitiesPages = GTK_WIDGET(gtk_builder_get_object(builder, "lbl_probabilitiesPages"));

  gtk_window_maximize(GTK_WINDOW(window));
  gtk_builder_connect_signals(builder, NULL);
  createProbabilitiesGrid();

  
  gtk_widget_show(window);
  gtk_main();
  return 0;
}

void initRCPairs(){
  initArray(rcArray, struct RCPair*, 16);
  for(int i = 0; i < rcArray.size; i++){
    rcArray.data[i] = NULL;
  }
}

gboolean on_da_chromosomes_draw(GtkWidget *widget, cairo_t *cr, gpointer user_data){
  do_drawing(cr);

  return FALSE;
}

void do_drawing(cairo_t *cr){
  double xCoordChromosome, yCoordGene;
  if(needsDrawing){
    needsDrawing = false;
    printf("SE LLAMA\n");
    if(areGenesInitialized()){
      Array_ints orders = getOrders();
      for(int i = 0; i < getTotalPossibleOrders(); i++){
        for(int k = 0; k < orders.data[i].used; k++){
          printf("%d,", orders.data[i].data[k]);
        }
        printf("\n");
        xCoordChromosome = 10+i*20;
        Array_double distances = getRelativeDistances(orders.data[i]);
        cairo_scale(cr, scale, scale);
        cairo_set_line_width(cr, 5);
        cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND); 
        cairo_move_to(cr, xCoordChromosome, 10);
        cairo_line_to(cr, xCoordChromosome, 110);
        cairo_move_to(cr, xCoordChromosome-10, 20);
        cairo_line_to(cr, xCoordChromosome+10, 20);
        for(int j = 0; j < distances.used; j++){
          printf("xcc:%lf, d:%lf,", xCoordChromosome, distances.data[j]);
          yCoordGene = 20+100*distances.data[j];
          cairo_move_to(cr, xCoordChromosome-10, yCoordGene);
          cairo_line_to(cr, xCoordChromosome+10, yCoordGene);
        }
        printf("\n");
        cairo_stroke(cr);
      }
    }
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
        gtk_widget_show(cell);
        if(row != 0 && col != 0){
          struct RCPair *rcp = malloc(sizeof(struct RCPair));
          rcp->row = row-1;
          rcp->col = col-1;
          appendArray(rcArray, struct RCPair*, rcp);
          g_signal_connect(cell, "insert-text", G_CALLBACK(on_txt_probability_insert_text), rcArray.data[rcArray.used-1]);
        }
      }
    }
  }
  gtk_grid_set_row_spacing(GTK_GRID(grid_probabilities), 10);
  gtk_grid_set_column_spacing(GTK_GRID(grid_probabilities), 10);
}

void fillProbabilitiesGrid(int total){
  if(totalProbs != total){
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
        } else {
          gtk_entry_set_text(GTK_ENTRY(cell), "");
        }
      } else if(col == 0 && row > 0){
        if(lRow-1 < getTotalGenes()){
          Array_char geneName = getGeneName(lRow-1);
          gtk_entry_set_text(GTK_ENTRY(cell), geneName.data);
        } else {
          gtk_entry_set_text(GTK_ENTRY(cell), "");
        }
      } else if(row > 0 && col > 0){
        if(lCol <= lRow){
          gtk_widget_set_sensitive(cell, FALSE);
        } else {
          gtk_widget_set_sensitive(cell, TRUE);
        }
        // Array_char d = getDescendant(probRow,probCol);
        // if(d.data == NULL){
        //   gtk_entry_set_text(GTK_ENTRY(cell), "");
        // } else {
        //   gtk_entry_set_text(GTK_ENTRY(cell), d.data);
        // }
      }
    }
  }
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
  for(int i = 0; i < rcArray.size; i++){
    if(rcArray.data[i]){
      free(rcArray.data[i]);
    }
  }
  freeArray(rcArray);
  if(areGenesInitialized()){
    freeGenes();
  }
  g_object_unref(builder);
  gtk_main_quit();
}

void on_btn_chooseFile_file_set(GtkFileChooserButton *f){
  gchar *fileName = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(f));
  printf("%s\n",fileName);
  fileChosen = true;
}

void on_btn_nFeatures_clicked(GtkButton *b){
  const gchar* val = gtk_entry_get_text(GTK_ENTRY(txt_nFeatures));
  int n = 0;
  if(gtk_entry_get_text_length(GTK_ENTRY(txt_nFeatures)) > 0){
    sscanf(val, "%d", &n);
  }
  if(n > 0){
    if(!areGenesInitialized()){
      initializeGenes(n);
    }
    createGeneNames(n);
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

void on_txt_probability_insert_text(GtkEditable* self, gchar* new_text, gint new_text_length, gint* position, gpointer user_data){
  double prob = 0.0;
  int row, col;
  if(new_text_length <= 1){
    if(*position == 0){
      if(new_text[0] != '0'){
        g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
      }
    } else if(*position == 1){
      if(new_text[0] != ','){
        g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
      }
    } else if(*position == 2){
      if(new_text[0] < '0' || new_text[0] >= '5'){
        g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
      } else {
        const gchar* c = gtk_entry_get_text(GTK_ENTRY(self));
        char result[18];
        strcpy(result,c);
        strcat(result,new_text);
        sscanf(result, "%lf", &prob);
        row = (*(struct RCPair*)user_data).row;
        col = (*(struct RCPair*)user_data).col;
        insertProbability(prob, row, col);
      }
    } else if(new_text[0] < '0' || new_text[0] > '9'){
      g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
    }
  } else {
      sscanf(new_text, "%lf", &prob);
      if(prob <= 0.0 || prob >= 0.5){
        g_signal_stop_emission_by_name (G_OBJECT(self), "insert_text");
      } else {
        row = (*(struct RCPair*)user_data).row;
        col = (*(struct RCPair*)user_data).col;
        insertProbability(prob, row, col);
      }
    }
}

void on_btn_zoomIn_clicked(GtkButton *b){
  scale += 0.1;
}

void on_btn_zoomOut_clicked(GtkButton *b){
  scale -= 0.1;
}

void on_btn_inferProbabilities_clicked(GtkButton *b){
  inferProbabilities();
}

void on_btn_draw_clicked(GtkButton *b){
  createOrders();
  needsDrawing = true;
  gtk_widget_queue_draw(window);
}