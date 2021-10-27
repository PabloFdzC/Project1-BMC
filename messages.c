#include "messages.h"

GtkBuilder *builderMessages;

void messagesWindow(GtkWidget *window, char m[]){
  builderMessages = gtk_builder_new();
  gtk_builder_add_from_file(builderMessages, "glade/messages.glade", NULL);

  GtkWidget *messages_window = GTK_WIDGET(gtk_builder_get_object(builderMessages, "messages_window"));
  gtk_window_resize (GTK_WINDOW(messages_window),350, 200);
  gtk_window_set_transient_for(GTK_WINDOW(messages_window), GTK_WINDOW(window));
  gtk_window_set_modal(GTK_WINDOW(messages_window), TRUE);
  gtk_builder_connect_signals(builderMessages, NULL);

  GtkWidget *lbl_message = GTK_WIDGET(gtk_builder_get_object(builderMessages, "lbl_message"));
  gtk_label_set_text(GTK_LABEL(lbl_message), m);  
  gtk_widget_show(messages_window);
  gtk_main();  
}

void on_messages_window_destroy(){
  g_object_unref(builderMessages);
  gtk_main_quit();
}