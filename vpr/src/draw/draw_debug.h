#pragma once
/**
 * @file draw_debug.h
 * 
 *  This file contains all functions regarding the graphics related to the setting of place and route breakpoints.
 * Manages creation of new Gtk Windows with debug options on use of the "Debug" button.
 */

#ifndef NO_GRAPHICS

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <gtk/gtk.h>
#include <string>
#include "breakpoint_state_globals.h"

/** debugger functions **/
void draw_debug_window();
void refresh_bpList();
void add_to_bpList(std::string bpDescription);
void set_moves_button_callback(GtkWidget* /*widget*/, GtkWidget* grid);
void set_temp_button_callback(GtkWidget* /*widgets*/, GtkWidget* grid);
void set_block_button_callback(GtkWidget* /*widget*/, GtkWidget* grid);
void set_router_iter_button_callback(GtkWidget* /*widget*/, GtkWidget* grid);
void set_net_id_button_callback(GtkWidget* /*widget*/, GtkWidget* grid);
void checkbox_callback(GtkWidget* widget);
void delete_bp_callback(GtkWidget* widget);
void advanced_button_callback();
void set_expression_button_callback(GtkWidget* /*widget*/, GtkWidget* grid);
void close_debug_window();
void close_advanced_window();
void ok_close_window(GtkWidget* /*widget*/, GtkWidget* window);
void invalid_breakpoint_entry_window(std::string error);
bool valid_expression(std::string exp);
void breakpoint_info_window(std::string bpDescription, BreakpointState draw_breakpoint_state, bool in_placer);

#endif /*NO_GRAPHICS*/
