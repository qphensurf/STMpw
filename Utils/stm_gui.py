#!/usr/bin/env python3

import os
import gzip
import json
import re
import time

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
except ImportError:
    tk = None
    filedialog = None
    messagebox = None
    ttk = None
    FigureCanvasTkAgg = None

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.widgets import Button, TextBox, RadioButtons
from matplotlib.transforms import Bbox
from scipy.interpolate import RegularGridInterpolator


class STMData:
    def __init__(self):
        self.filename = None
        self.N_sampling_z = None
        self.ngy = None
        self.ngx = None
        self.z = None
        self.y = None
        self.x = None
        self.intensity = None

    def read(self, filename):
        self.filename = filename
        open_func = gzip.open if filename.endswith(".gz") else open

        with open_func(filename, "rt") as f:
            first = f.readline().split()
            if len(first) < 3:
                raise ValueError("First line must contain: N_sampling_z ngy ngx")
            self.N_sampling_z, self.ngy, self.ngx = map(int, first[:3])
            self.z = np.zeros(self.N_sampling_z, dtype=float)
            self.y = np.zeros(self.ngy, dtype=float)
            self.x = np.zeros(self.ngx, dtype=float)
            self.intensity = np.zeros((self.ngx, self.ngy, self.N_sampling_z), dtype=float)
            expected = self.N_sampling_z * self.ngy * self.ngx
            nread = 0
            for iz in range(self.N_sampling_z):
                for iy in range(self.ngy):
                    for ix in range(self.ngx):
                        line = f.readline()
                        if not line:
                            raise ValueError(f"Unexpected end of file after {nread} data rows. Expected {expected}.")
                        parts = line.split()
                        if len(parts) < 4:
                            raise ValueError(f"Malformed data row after {nread} rows: {line}")
                        zz = float(parts[0])
                        yy = float(parts[1])
                        xx = float(parts[2])
                        val = float(parts[3])
                        self.z[iz] = zz
                        self.y[iy] = yy
                        self.x[ix] = xx
                        self.intensity[ix, iy, iz] = abs(val)
                        nread += 1
        self._check_axes()

    def _check_axes(self):
        if len(np.unique(self.x)) != len(self.x):
            raise ValueError("Repeated values detected in x grid.")
        if len(np.unique(self.y)) != len(self.y):
            raise ValueError("Repeated values detected in y grid.")
        if len(np.unique(self.z)) != len(self.z):
            raise ValueError("Repeated values detected in z grid.")
        if not np.all(np.diff(self.x) > 0):
            order = np.argsort(self.x)
            self.x = self.x[order]
            self.intensity = self.intensity[order, :, :]
        if not np.all(np.diff(self.y) > 0):
            order = np.argsort(self.y)
            self.y = self.y[order]
            self.intensity = self.intensity[:, order, :]
        if not np.all(np.diff(self.z) > 0):
            order = np.argsort(self.z)
            self.z = self.z[order]
            self.intensity = self.intensity[:, :, order]


def constant_current_map_linear(intensity, z, current):
    ngx, ngy, nz = intensity.shape
    zmap = np.full((ngx, ngy), np.nan, dtype=float)
    for ix in range(ngx):
        for iy in range(ngy):
            prof = intensity[ix, iy, :]
            for iz in range(nz - 1):
                i0 = prof[iz]
                i1 = prof[iz + 1]
                if not np.isfinite(i0) or not np.isfinite(i1) or i0 == i1:
                    continue
                if (current - i0) * (current - i1) <= 0.0:
                    z0 = z[iz]
                    z1 = z[iz + 1]
                    zmap[ix, iy] = z0 + (current - i0) * (z1 - z0) / (i1 - i0)
                    break
    return zmap


def constant_current_map_log(intensity, z, current):
    if current <= 0.0:
        raise ValueError("Log interpolation requires a positive current setpoint.")
    ngx, ngy, nz = intensity.shape
    zmap = np.full((ngx, ngy), np.nan, dtype=float)
    log_current = np.log(current)
    for ix in range(ngx):
        for iy in range(ngy):
            prof = intensity[ix, iy, :]
            for iz in range(nz - 1):
                i0 = prof[iz]
                i1 = prof[iz + 1]
                if i0 <= 0.0 or i1 <= 0.0 or not np.isfinite(i0) or not np.isfinite(i1) or i0 == i1:
                    continue
                if (current - i0) * (current - i1) <= 0.0:
                    li0 = np.log(i0)
                    li1 = np.log(i1)
                    if li0 == li1:
                        continue
                    z0 = z[iz]
                    z1 = z[iz + 1]
                    zmap[ix, iy] = z0 + (log_current - li0) * (z1 - z0) / (li1 - li0)
                    break
    return zmap


def constant_height_map(intensity, z, height):
    if height < np.min(z) or height > np.max(z):
        raise ValueError(f"Requested height {height} is outside z range [{np.min(z)}, {np.max(z)}].")
    ngx, ngy, nz = intensity.shape
    imap = np.full((ngx, ngy), np.nan, dtype=float)
    for ix in range(ngx):
        for iy in range(ngy):
            imap[ix, iy] = np.interp(height, z, intensity[ix, iy, :])
    return imap


def periodic_step(axis):
    if len(axis) < 2:
        raise ValueError("Axis must contain at least two points.")
    return float(np.mean(np.diff(axis)))


def extend_periodic_image(image, x, y, nx, ny):
    if nx < 1 or ny < 1:
        raise ValueError("nx and ny must be positive integers.")
    dx = periodic_step(x)
    dy = periodic_step(y)
    lx = (x[-1] - x[0]) + dx
    ly = (y[-1] - y[0]) + dy
    extended_image = np.tile(image, (nx, ny))
    extended_x = np.concatenate([x + ix * lx for ix in range(nx)])
    extended_y = np.concatenate([y + iy * ly for iy in range(ny)])
    return extended_image, extended_x, extended_y, lx, ly


def axis_edges_from_centers(axis):
    axis = np.asarray(axis, dtype=float)
    if len(axis) < 2:
        raise ValueError("Axis must contain at least two points.")
    edges = np.empty(len(axis) + 1, dtype=float)
    edges[1:-1] = 0.5 * (axis[:-1] + axis[1:])
    edges[0] = axis[0] - 0.5 * (axis[1] - axis[0])
    edges[-1] = axis[-1] + 0.5 * (axis[-1] - axis[-2])
    return edges


def normalized_vector(xval, yval, name):
    vec = np.array([float(xval), float(yval)], dtype=float)
    norm = np.linalg.norm(vec)
    if norm == 0.0:
        raise ValueError(f"{name} vector cannot have zero length.")
    return vec / norm


def make_lattice_transform(a1_x, a1_y, a2_x, a2_y):
    ex = normalized_vector(a1_x, a1_y, "a1")
    ey = normalized_vector(a2_x, a2_y, "a2")
    basis = np.array([[ex[0], ey[0]], [ex[1], ey[1]]], dtype=float)
    det = np.linalg.det(basis)
    if abs(det) < 1.0e-12:
        raise ValueError("a1 and a2 directions are linearly dependent.")
    return basis, np.linalg.inv(basis)


def lattice_to_cartesian(u, v, basis):
    return basis[0, 0] * u + basis[0, 1] * v, basis[1, 0] * u + basis[1, 1] * v


def cartesian_to_lattice(xcart, ycart, inverse_basis):
    return inverse_basis[0, 0] * xcart + inverse_basis[0, 1] * ycart, inverse_basis[1, 0] * xcart + inverse_basis[1, 1] * ycart


def remove_nonfinite_profile_points(distance, profile):
    distance = np.asarray(distance, dtype=float)
    profile = np.asarray(profile, dtype=float)
    mask = np.isfinite(distance) & np.isfinite(profile)
    return distance[mask], profile[mask]


def point_to_segment_distance(px, py, x0, y0, x1, y1):
    vx = x1 - x0
    vy = y1 - y0
    wx = px - x0
    wy = py - y0
    c1 = vx * wx + vy * wy
    if c1 <= 0.0:
        return np.sqrt((px - x0) ** 2 + (py - y0) ** 2)
    c2 = vx * vx + vy * vy
    if c2 <= 0.0:
        return np.sqrt((px - x0) ** 2 + (py - y0) ** 2)
    if c1 >= c2:
        return np.sqrt((px - x1) ** 2 + (py - y1) ** 2)
    b = c1 / c2
    bx = x0 + b * vx
    by = y0 + b * vy
    return np.sqrt((px - bx) ** 2 + (py - by) ** 2)


def load_gnuplot_palette_file(filename, cmap_name="custom_gnuplot_palette"):
    with open(filename, "r") as f:
        text = f.read()
    text = text.replace("\\", " ").replace("\n", " ")
    match = re.search(r"defined\s*\((.*)\)", text, re.IGNORECASE)
    if match is None:
        raise ValueError("Could not find a 'set palette defined(...)' block in the palette file.")
    points = []
    for entry in match.group(1).split(","):
        parts = entry.split()
        if len(parts) < 4:
            continue
        try:
            xval = float(parts[0])
            rval = min(max(float(parts[1]), 0.0), 1.0)
            gval = min(max(float(parts[2]), 0.0), 1.0)
            bval = min(max(float(parts[3]), 0.0), 1.0)
        except ValueError:
            continue
        points.append((xval, rval, gval, bval))
    if len(points) < 2:
        raise ValueError("The palette file must contain at least two valid color points.")
    points.sort(key=lambda item: item[0])
    xmin = points[0][0]
    xmax = points[-1][0]
    if xmax == xmin:
        raise ValueError("Palette coordinate range is zero.")
    color_list = [((xval - xmin) / (xmax - xmin), (rval, gval, bval)) for xval, rval, gval, bval in points]
    return LinearSegmentedColormap.from_list(cmap_name, color_list)


class CustomPaletteDialog:
    def __init__(self, parent_gui, initial_directory=".", previous_palette_name="gray", previous_cmap="gray"):
        self.parent_gui = parent_gui
        self.previous_palette_name = previous_palette_name
        self.previous_cmap = previous_cmap
        self.selected_palette_file = None
        self.selected_cmap = None
        self.selected_palette_name = None
        if tk._default_root is None:
            self.root = tk.Tk()
        else:
            self.root = tk.Toplevel(tk._default_root)
        self.root.title("Load custom palette from file")
        self.root.geometry("850x520")
        self.current_directory = tk.StringVar(value=os.path.abspath(initial_directory))
        self.status_text = tk.StringVar(value="Select a palette file to preview it.")
        self._build_dialog()
        self._populate_file_list()
        self.root.protocol("WM_DELETE_WINDOW", self._cancel)

    def _build_dialog(self):
        top_frame = tk.Frame(self.root)
        top_frame.pack(side=tk.TOP, fill=tk.X, padx=8, pady=8)
        tk.Label(top_frame, text="Directory:").pack(side=tk.LEFT)
        self.directory_entry = tk.Entry(top_frame, textvariable=self.current_directory, width=70)
        self.directory_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        tk.Button(top_frame, text="Browse", command=self._browse_directory).pack(side=tk.LEFT, padx=5)
        tk.Button(top_frame, text="Refresh", command=self._populate_file_list).pack(side=tk.LEFT)
        main_frame = tk.Frame(self.root)
        main_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=8, pady=8)
        left_frame = tk.Frame(main_frame)
        left_frame.pack(side=tk.LEFT, fill=tk.Y)
        tk.Label(left_frame, text="Palette files").pack(anchor=tk.W)
        self.file_listbox = tk.Listbox(left_frame, width=34, height=20)
        self.file_listbox.pack(side=tk.LEFT, fill=tk.Y)
        scrollbar = tk.Scrollbar(left_frame, orient=tk.VERTICAL)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.file_listbox.config(yscrollcommand=scrollbar.set)
        scrollbar.config(command=self.file_listbox.yview)
        self.file_listbox.bind("<<ListboxSelect>>", self._on_file_selected)
        right_frame = tk.Frame(main_frame)
        right_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=15)
        tk.Label(right_frame, text="Palette preview").pack(anchor=tk.W)
        self.preview_figure, self.preview_ax = plt.subplots(figsize=(6.0, 1.6))
        self.preview_canvas = FigureCanvasTkAgg(self.preview_figure, master=right_frame)
        self.preview_canvas.get_tk_widget().pack(fill=tk.X, pady=10)
        self.preview_ax.set_yticks([])
        self.preview_ax.set_xlabel("Normalized palette coordinate")
        self.preview_ax.set_title("No palette selected")
        self.preview_canvas.draw()
        tk.Label(right_frame, textvariable=self.status_text, anchor=tk.W, justify=tk.LEFT).pack(fill=tk.X, pady=10)
        button_frame = tk.Frame(self.root)
        button_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=8, pady=8)
        tk.Button(button_frame, text="Cancel", command=self._cancel, width=12).pack(side=tk.RIGHT, padx=5)
        tk.Button(button_frame, text="Apply", command=self._apply, width=12).pack(side=tk.RIGHT, padx=5)

    def _browse_directory(self):
        directory = filedialog.askdirectory(initialdir=self.current_directory.get(), title="Select palette directory")
        if directory:
            self.current_directory.set(directory)
            self._populate_file_list()

    def _populate_file_list(self):
        self.file_listbox.delete(0, tk.END)
        directory = self.current_directory.get()
        if not os.path.isdir(directory):
            self.status_text.set(f"Directory does not exist: {directory}")
            return
        candidates = []
        for filename in sorted(os.listdir(directory)):
            full_path = os.path.join(directory, filename)
            if not os.path.isfile(full_path):
                continue
            lower_name = filename.lower()
            if lower_name.endswith((".pal", ".plt", ".gp", ".gnu", ".txt")):
                candidates.append(filename)
        for filename in candidates:
            self.file_listbox.insert(tk.END, filename)
        self.status_text.set("Click a palette file to preview and temporarily apply it." if candidates else "No palette-like files found in this directory.")

    def _on_file_selected(self, event=None):
        selection = self.file_listbox.curselection()
        if not selection:
            return
        filename = self.file_listbox.get(selection[0])
        full_path = os.path.join(self.current_directory.get(), filename)
        try:
            cmap_name = f"custom_{os.path.splitext(os.path.basename(full_path))[0]}"
            cmap = load_gnuplot_palette_file(full_path, cmap_name=cmap_name)
            self.selected_palette_file = full_path
            self.selected_cmap = cmap
            self.selected_palette_name = f"custom:{full_path}"
            self._draw_preview(cmap, full_path)
            self.parent_gui._apply_palette_temporarily(self.selected_palette_name, self.selected_cmap)
            self.status_text.set(f"Previewing and temporarily applying: {full_path}")
        except Exception as exc:
            self.status_text.set(f"Error loading palette: {exc}")
            if messagebox is not None:
                messagebox.showerror("Palette error", str(exc))

    def _draw_preview(self, cmap, filename):
        self.preview_ax.clear()
        gradient = np.linspace(0.0, 1.0, 512)
        gradient = np.vstack((gradient, gradient))
        self.preview_ax.imshow(gradient, aspect="auto", cmap=cmap, origin="lower", extent=[0.0, 1.0, 0.0, 1.0])
        self.preview_ax.set_yticks([])
        self.preview_ax.set_xlabel("Normalized palette coordinate")
        self.preview_ax.set_title(os.path.basename(filename))
        self.preview_canvas.draw()

    def _apply(self):
        if self.selected_cmap is None or self.selected_palette_name is None:
            if messagebox is not None:
                messagebox.showwarning("No palette selected", "Select a palette file before pressing Apply.")
            return
        self.parent_gui._finalize_custom_palette(self.selected_palette_name, self.selected_palette_file, self.selected_cmap)
        self.root.destroy()

    def _cancel(self):
        self.parent_gui._restore_palette_after_cancel(self.previous_palette_name, self.previous_cmap)
        self.root.destroy()

    def show(self):
        # Do not use mainloop() here. If a hidden Tk root was already created
        # for clipboard paste, mainloop() may keep running after this Toplevel
        # is destroyed. wait_window() returns when this dialog is closed.
        try:
            self.root.grab_set()
        except Exception:
            pass
        self.root.wait_window(self.root)


class STMLineProfileGUI:
    def __init__(self):
        self.data = STMData()
        self.mode = "constant current"
        self.interpolation = "log"
        self.current_value = 1.0e-10
        self.height_value = 10.0
        self.nx = 1
        self.ny = 1
        self.a1_x = 1.0
        self.a1_y = 0.0
        self.a2_x = 0.0
        self.a2_y = 1.0
        self.basis = np.eye(2)
        self.inverse_basis = np.eye(2)
        self.cart_x_edges = None
        self.cart_y_edges = None
        self.palette = "gray"
        self.active_cmap = "gray"
        self.custom_palette_file = None
        self.custom_palette_cmap = None
        self.palette_button = None
        self._ignore_palette_callback = False
        self.palette_options = ("gray", "viridis", "plasma", "inferno", "magma", "cividis", "hot", "afmhot", "turbo", "seismic", "load from file...")
        self.base_image = None
        self.image = None
        self.extended_x = None
        self.extended_y = None
        self.cell_lx = None
        self.cell_ly = None
        self.image_label = ""
        self.image_units = ""
        self.image_interpolator = None
        self.lines = []
        self.profile_counter = 0
        self.selected_label = None
        self.pick_start = None
        self.preview_line = None
        self.pending_action = None
        self.last_mouse_xy = None
        self.drag_action = None
        self.drag_line_label = None
        self.drag_start_mouse = None
        self.drag_original_line = None
        self.live_profile_line = None
        self.live_profile_label = "preview"
        self.last_live_profile_update = 0.0
        self.live_profile_update_interval = 0.05
        self.fig_img = None
        self.ax_img = None
        self.ax_cbar = None
        self.cbar_position = [0.635, 0.12, 0.018, 0.80]
        self.im = None
        self.cbar = None
        self.coord_text = None
        self.status_text = None
        self.fig_prof = None
        self.ax_prof = None
        self.textbox_tab_order = []
        self.active_textbox = None
        self._clipboard_root = None
        self.session_button = None
        self.save_clean_image_button = None
        self.undo_button = None
        self.redo_button = None
        self.undo_stack = []
        self.redo_stack = []
        self.undo_limit = 50
        self._restoring_history = False
        self._build_gui()

    def _build_gui(self):
        self.fig_img = plt.figure(figsize=(13.4, 9.4))
        self.fig_img.canvas.manager.set_window_title("STM image and line selection")
        self.ax_img = self.fig_img.add_axes([0.06, 0.12, 0.56, 0.80])
        self.ax_cbar = self.fig_img.add_axes(self.cbar_position)
        self.ax_cbar.set_visible(False)
        self.ax_img.set_title("Load STMpw data file")
        self.ax_img.set_xlabel("X (Å)")
        self.ax_img.set_ylabel("Y (Å)")
        self.ax_img.set_aspect("equal")
        ax_file = self.fig_img.add_axes([0.755, 0.925, 0.205, 0.035])
        self.file_box = TextBox(ax_file, "File", initial="TH.dat")
        ax_mode = self.fig_img.add_axes([0.755, 0.810, 0.205, 0.085])
        self.mode_radio = RadioButtons(ax_mode, ("constant current", "constant height"))
        ax_interp = self.fig_img.add_axes([0.755, 0.715, 0.205, 0.065])
        self.interp_radio = RadioButtons(ax_interp, ("log", "linear"))
        ax_current = self.fig_img.add_axes([0.755, 0.655, 0.205, 0.035])
        self.current_box = TextBox(ax_current, "Current", initial="1e-10")
        ax_height = self.fig_img.add_axes([0.755, 0.605, 0.205, 0.035])
        self.height_box = TextBox(ax_height, "Height Å", initial="10")
        ax_nx = self.fig_img.add_axes([0.755, 0.555, 0.090, 0.035])
        self.nx_box = TextBox(ax_nx, "nx", initial="1")
        ax_ny = self.fig_img.add_axes([0.870, 0.555, 0.090, 0.035])
        self.ny_box = TextBox(ax_ny, "ny", initial="1")
        ax_a1x = self.fig_img.add_axes([0.755, 0.510, 0.090, 0.032])
        self.a1x_box = TextBox(ax_a1x, "a1_x", initial="1.0")
        ax_a1y = self.fig_img.add_axes([0.870, 0.510, 0.090, 0.032])
        self.a1y_box = TextBox(ax_a1y, "a1_y", initial="0.0")
        ax_a2x = self.fig_img.add_axes([0.755, 0.472, 0.090, 0.032])
        self.a2x_box = TextBox(ax_a2x, "a2_x", initial="0.0")
        ax_a2y = self.fig_img.add_axes([0.870, 0.472, 0.090, 0.032])
        self.a2y_box = TextBox(ax_a2y, "a2_y", initial="1.0")
        self.fig_img.text(0.755, 0.452, "Palette", fontsize=10)
        ax_palette = self.fig_img.add_axes([0.755, 0.410, 0.205, 0.040])
        self.palette_button = Button(ax_palette, self._palette_button_label())
        ax_load = self.fig_img.add_axes([0.755, 0.350, 0.205, 0.035])
        self.load_button = Button(ax_load, "Load/plot image")
        ax_line_id = self.fig_img.add_axes([0.755, 0.305, 0.085, 0.035])
        self.line_id_box = TextBox(ax_line_id, "Line", initial="")
        ax_select_line = self.fig_img.add_axes([0.860, 0.305, 0.100, 0.035])
        self.select_line_button = Button(ax_select_line, "Select line")
        ax_pick_nearest = self.fig_img.add_axes([0.755, 0.260, 0.205, 0.035])
        self.pick_nearest_button = Button(ax_pick_nearest, "Pick nearest line")
        ax_modify = self.fig_img.add_axes([0.755, 0.215, 0.100, 0.035])
        self.modify_button = Button(ax_modify, "Modify selected")
        ax_delete = self.fig_img.add_axes([0.860, 0.215, 0.100, 0.035])
        self.delete_button = Button(ax_delete, "Delete selected")
        ax_clear = self.fig_img.add_axes([0.755, 0.170, 0.100, 0.035])
        self.clear_button = Button(ax_clear, "Clear lines")
        ax_save_data = self.fig_img.add_axes([0.860, 0.170, 0.100, 0.035])
        self.save_data_button = Button(ax_save_data, "Save profiles")
        ax_save_fig = self.fig_img.add_axes([0.755, 0.125, 0.100, 0.035])
        self.save_fig_button = Button(ax_save_fig, "Save image")

        ax_save_clean = self.fig_img.add_axes([0.860, 0.125, 0.100, 0.035])
        self.save_clean_image_button = Button(ax_save_clean, "Save clean")

        ax_undo = self.fig_img.add_axes([0.755, 0.080, 0.100, 0.035])
        self.undo_button = Button(ax_undo, "Undo")

        ax_redo = self.fig_img.add_axes([0.860, 0.080, 0.100, 0.035])
        self.redo_button = Button(ax_redo, "Redo")

        ax_session = self.fig_img.add_axes([0.755, 0.035, 0.070, 0.035])
        self.session_button = Button(ax_session, "Session")

        ax_help = self.fig_img.add_axes([0.835, 0.035, 0.055, 0.035])
        self.help_button = Button(ax_help, "Help")

        ax_close = self.fig_img.add_axes([0.900, 0.035, 0.060, 0.035])
        self.close_button = Button(ax_close, "Close")
        self.coord_text = self.fig_img.text(0.06, 0.050, "Mouse: X = -, Y = -, value = -", fontsize=9)
        self.status_text = self.fig_img.text(0.06, 0.025, "Status: draw mode. Click first point, then second point.", fontsize=9)
        self.mode_radio.on_clicked(self._set_mode)
        self.interp_radio.on_clicked(self._set_interpolation)
        self.palette_button.on_clicked(self._open_palette_dropdown_dialog)
        self.load_button.on_clicked(self._load_and_plot)
        self.select_line_button.on_clicked(self._select_line_from_textbox)
        self.pick_nearest_button.on_clicked(self._start_pick_nearest_line)
        self.modify_button.on_clicked(self._start_modify_selected_line)
        self.delete_button.on_clicked(self._delete_selected_line)
        self.clear_button.on_clicked(self._clear_lines)
        self.save_data_button.on_clicked(self._save_profiles)
        self.save_fig_button.on_clicked(self._save_figures)
        self.save_clean_image_button.on_clicked(self._save_clean_image)
        self.session_button.on_clicked(self._open_session_dialog)
        self.undo_button.on_clicked(self._undo)
        self.redo_button.on_clicked(self._redo)
        self.help_button.on_clicked(self._show_help)
        self.close_button.on_clicked(self._close_program)
        self._register_textboxes_for_editing()
        self.fig_img.canvas.mpl_connect("key_press_event", self._session_key_press)
        self.fig_img.canvas.mpl_connect("key_press_event", self._textbox_key_press)
        self.fig_img.canvas.mpl_connect("button_press_event", self._textbox_mouse_press)
        self.fig_img.canvas.mpl_connect("button_press_event", self._on_click)
        self.fig_img.canvas.mpl_connect("button_release_event", self._on_release)
        self.fig_img.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self._ensure_profile_window()

    def _format_status_for_display(self, text, max_chars=155):
        text = str(text).replace("\n", " ")
        if len(text) <= max_chars:
            return text

        head = max_chars // 2 - 4
        tail = max_chars - head - 5
        shortened = text[:head] + " ... " + text[-tail:]
        print(f"Full status: {text}")
        return shortened

    def _set_status(self, text):
        compact_text = self._format_status_for_display(text)
        self.status_text.set_text(f"Status: {compact_text}")
        self.fig_img.canvas.draw_idle()

    def _on_profile_window_closed(self, event=None):
        self.live_profile_line = None
        self.fig_prof = None
        self.ax_prof = None

    def _profile_window_is_available(self):
        if self.fig_prof is None or self.ax_prof is None:
            return False

        try:
            figure_number = self.fig_prof.number
        except Exception:
            return False

        if not plt.fignum_exists(figure_number):
            return False

        try:
            manager = self.fig_prof.canvas.manager
        except Exception:
            return False

        return manager is not None

    def _show_profile_window_now(self):
        if self.fig_prof is None:
            return

        try:
            self.fig_prof.show()
        except Exception:
            pass

        try:
            manager = self.fig_prof.canvas.manager
            if manager is not None and hasattr(manager, "show"):
                manager.show()
        except Exception:
            pass

        try:
            self.fig_prof.canvas.draw_idle()
            self.fig_prof.canvas.flush_events()
        except Exception:
            pass

    def _ensure_profile_window(self):
        if self._profile_window_is_available():
            self._show_profile_window_now()
            return

        self.fig_prof, self.ax_prof = plt.subplots(figsize=(8, 5))
        self.fig_prof.canvas.manager.set_window_title("Line profiles")
        self.fig_prof.canvas.mpl_connect("close_event", self._on_profile_window_closed)
        self.ax_prof.set_xlabel("Distance (Å)")
        self.ax_prof.set_ylabel("Profile value")
        self.ax_prof.set_title("Line profiles")
        self.ax_prof.grid(True)
        self._show_profile_window_now()

    def _clear_live_profile_preview(self, redraw=True):
        if self.live_profile_line is not None:
            try:
                self.live_profile_line.remove()
            except Exception:
                pass
            self.live_profile_line = None

        if redraw and self._profile_window_is_available():
            try:
                self.fig_prof.canvas.draw_idle()
            except Exception:
                pass

    def _update_live_profile_preview(self, x0, y0, x1, y1, force=False):
        now = time.monotonic()
        if not force and (now - self.last_live_profile_update) < self.live_profile_update_interval:
            return
        self.last_live_profile_update = now

        if self.image is None or self.image_interpolator is None:
            return

        self._ensure_profile_window()

        try:
            distance, profile = self._extract_profile(x0, y0, x1, y1)
        except Exception:
            return

        if len(distance) < 2:
            return

        if self.live_profile_line is None or self.live_profile_line.axes is not self.ax_prof:
            line_obj, = self.ax_prof.plot(
                distance,
                profile,
                linestyle="--",
                linewidth=2.0,
                color="red",
                label=self.live_profile_label
            )
            self.live_profile_line = line_obj
        else:
            self.live_profile_line.set_data(distance, profile)

        try:
            self.ax_prof.relim()
            self.ax_prof.autoscale_view()
            legend = self.ax_prof.get_legend()
            if legend is None:
                self.ax_prof.legend()
            self.fig_prof.canvas.draw_idle()
        except Exception:
            pass

    def _line_length_angle_text(self, x0, y0, x1, y1):
        dx = float(x1 - x0)
        dy = float(y1 - y0)
        length = float(np.sqrt(dx * dx + dy * dy))
        angle = float(np.degrees(np.arctan2(dy, dx)))
        return f"length = {length:.6g} Å, angle = {angle:.3f}°"

    def _drag_hit_tolerance(self):
        if self.cart_x_edges is None or self.cart_y_edges is None:
            return 0.1

        xmin = float(np.nanmin(self.cart_x_edges))
        xmax = float(np.nanmax(self.cart_x_edges))
        ymin = float(np.nanmin(self.cart_y_edges))
        ymax = float(np.nanmax(self.cart_y_edges))
        diag = float(np.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2))

        if self.extended_x is not None and self.extended_y is not None and len(self.extended_x) > 1 and len(self.extended_y) > 1:
            step_x = abs(float(np.mean(np.diff(self.extended_x))))
            step_y = abs(float(np.mean(np.diff(self.extended_y))))
            grid_tol = 3.0 * min(step_x, step_y)
        else:
            grid_tol = 0.0

        return max(grid_tol, 0.015 * diag)

    def _copy_line_for_drag(self, line):
        copied = {}
        for key, value in line.items():
            if isinstance(value, np.ndarray):
                copied[key] = value.copy()
            else:
                copied[key] = value
        return copied

    def _restore_line_from_drag_copy(self):
        if self.drag_line_label is None or self.drag_original_line is None:
            return
        line = self._get_line_by_label(self.drag_line_label)
        if line is not None:
            line.clear()
            line.update(self._copy_line_for_drag(self.drag_original_line))

    def _clear_drag_state(self):
        self.drag_action = None
        self.drag_line_label = None
        self.drag_start_mouse = None
        self.drag_original_line = None

    def _hit_test_selected_line_for_drag(self, x, y):
        if self.selected_label is None:
            return None

        line = self._get_line_by_label(self.selected_label)
        if line is None:
            return None

        tolerance = self._drag_hit_tolerance()
        d0 = float(np.sqrt((x - line["x0"]) ** 2 + (y - line["y0"]) ** 2))
        d1 = float(np.sqrt((x - line["x1"]) ** 2 + (y - line["y1"]) ** 2))
        segment_distance = point_to_segment_distance(x, y, line["x0"], line["y0"], line["x1"], line["y1"])

        if d0 <= tolerance and d0 <= d1:
            return "start"
        if d1 <= tolerance:
            return "end"
        if segment_distance <= tolerance:
            return "line"

        return None

    def _start_line_drag(self, action, x, y):
        line = self._get_line_by_label(self.selected_label)
        if line is None:
            self._set_status("selected line not found")
            return False

        self._push_undo_state(f"drag modify {self.selected_label}")
        self.drag_action = action
        self.drag_line_label = self.selected_label
        self.drag_start_mouse = (x, y)
        self.drag_original_line = self._copy_line_for_drag(line)

        if action == "start":
            self._set_status(f"dragging first point of {self.selected_label}")
        elif action == "end":
            self._set_status(f"dragging second point of {self.selected_label}")
        else:
            self._set_status(f"dragging whole {self.selected_label}")

        return True

    def _update_dragged_line_coordinates(self, x, y):
        if self.drag_action is None or self.drag_original_line is None or self.drag_line_label is None:
            return

        line = self._get_line_by_label(self.drag_line_label)
        if line is None:
            return

        original = self.drag_original_line

        if self.drag_action == "start":
            line["x0"] = x
            line["y0"] = y
            line["x1"] = original["x1"]
            line["y1"] = original["y1"]
        elif self.drag_action == "end":
            line["x0"] = original["x0"]
            line["y0"] = original["y0"]
            line["x1"] = x
            line["y1"] = y
        elif self.drag_action == "line":
            start_x, start_y = self.drag_start_mouse
            dx = x - start_x
            dy = y - start_y
            new_x0, new_y0 = self._clamp_to_valid_domain(original["x0"] + dx, original["y0"] + dy)
            new_x1, new_y1 = self._clamp_to_valid_domain(original["x1"] + dx, original["y1"] + dy)
            line["x0"] = new_x0
            line["y0"] = new_y0
            line["x1"] = new_x1
            line["y1"] = new_y1

    def _finalize_dragged_line(self):
        if self.drag_line_label is None:
            return

        line = self._get_line_by_label(self.drag_line_label)
        if line is None:
            self._clear_drag_state()
            return

        distance, profile = self._extract_profile(line["x0"], line["y0"], line["x1"], line["y1"])
        if len(distance) < 2:
            print("Dragged line has fewer than two valid finite profile points. Drag was cancelled.")
            self._restore_line_from_drag_copy()
            self._clear_live_profile_preview(redraw=False)
            self._clear_drag_state()
            self._refresh_all_line_display()
            self._set_status("drag cancelled because line has too few valid points")
            return

        self._clear_live_profile_preview(redraw=False)
        line.update({
            "distance": distance,
            "profile": profile,
            "mode": self.mode,
            "current": self.current_value,
            "height": self.height_value,
            "units": self.image_units,
            "image_label": self.image_label,
            "nx": self.nx,
            "ny": self.ny,
            "palette": self.palette
        })

        label = self.drag_line_label
        self._clear_drag_state()
        self._refresh_all_line_display()
        self._show_selected_profile_statistics(f"modified {label}")

    def _register_textboxes_for_editing(self):
        names = ["file_box", "current_box", "height_box", "nx_box", "ny_box", "a1x_box", "a1y_box", "a2x_box", "a2y_box", "line_id_box"]
        self.textbox_tab_order = []
        for name in names:
            box = getattr(self, name, None)
            if box is not None:
                self.textbox_tab_order.append(box)

    def _textbox_mouse_press(self, event):
        for box in self.textbox_tab_order:
            if event.inaxes == box.ax:
                self.active_textbox = box
                if event.button == 2:
                    self._paste_into_textbox(box, prefer_primary=True)
                return

    def _textbox_key_press(self, event):
        if self.active_textbox is None:
            return
        key = str(event.key).lower() if event.key is not None else ""
        if key == "tab":
            self._focus_next_textbox(step=1)
            return
        if key in ("shift+tab", "shift tab"):
            self._focus_next_textbox(step=-1)
            return
        if key in ("ctrl+v", "control+v", "cmd+v", "super+v", "shift+insert", "ctrl+insert", "control+insert"):
            self._paste_into_textbox(self.active_textbox, prefer_primary=False)
            return

    def _focus_next_textbox(self, step=1):
        if not self.textbox_tab_order:
            return
        if self.active_textbox not in self.textbox_tab_order:
            new_box = self.textbox_tab_order[0]
        else:
            index = self.textbox_tab_order.index(self.active_textbox)
            new_box = self.textbox_tab_order[(index + step) % len(self.textbox_tab_order)]
        if self.active_textbox is not None:
            try:
                self.active_textbox.stop_typing()
            except Exception:
                pass
        self.active_textbox = new_box
        try:
            new_box.begin_typing()
            new_box.cursor_index = len(new_box.text)
            new_box._rendercursor()
        except Exception:
            pass
        self.fig_img.canvas.draw_idle()

    def _paste_into_textbox(self, box, prefer_primary=False):
        paste_text = self._get_clipboard_text(prefer_primary=prefer_primary)
        if paste_text is None:
            return
        paste_text = self._normalize_pasted_text(paste_text)
        if paste_text == "":
            return
        try:
            if not box.capturekeystrokes:
                box.begin_typing()
            cursor = getattr(box, "cursor_index", len(box.text))
            cursor = max(0, min(cursor, len(box.text)))
            new_text = box.text[:cursor] + paste_text + box.text[cursor:]
            box.cursor_index = cursor + len(paste_text)
            box.set_val(new_text)
            box._rendercursor()
            self.active_textbox = box
            self.fig_img.canvas.draw_idle()
        except Exception as exc:
            print(f"Could not paste into text box: {exc}")
            self._set_status(f"could not paste into text box: {exc}")

    def _get_clipboard_text(self, prefer_primary=False):
        if tk is None:
            print("tkinter is not available, clipboard paste is disabled.")
            self._set_status("tkinter is not available, clipboard paste is disabled")
            return None

        root = None
        try:
            root = tk.Tk()
            root.withdraw()

            if prefer_primary:
                try:
                    text = root.selection_get(selection="PRIMARY")
                    return self._normalize_pasted_text(text)
                except Exception:
                    pass

            try:
                text = root.clipboard_get()
                return self._normalize_pasted_text(text)
            except Exception:
                return None

        except Exception:
            return None
        finally:
            if root is not None:
                try:
                    root.destroy()
                except Exception:
                    pass

    def _normalize_pasted_text(self, text):
        text = str(text)
        text = text.replace("\r", " ").replace("\n", " ").replace("\t", " ")
        text = text.replace("\u00a0", " ")
        text = text.replace("\u2212", "-")
        text = text.replace("\u2013", "-")
        text = text.replace("\u2014", "-")
        return " ".join(text.split())

    def _close_program(self, event=None):
        """Close the STM GUI, profile plot, help window, palette dialog, and all Matplotlib windows."""
        try:
            if self.status_text is not None:
                self.status_text.set_text("Status: closing program")
                self.fig_img.canvas.draw_idle()
        except Exception:
            pass

        if tk is not None:
            root = getattr(tk, "_default_root", None)
            if root is not None:
                try:
                    for child in list(root.winfo_children()):
                        try:
                            child.destroy()
                        except Exception:
                            pass
                except Exception:
                    pass
                try:
                    root.destroy()
                except Exception:
                    pass
                try:
                    tk._default_root = None
                except Exception:
                    pass

        try:
            plt.close("all")
        except Exception:
            pass

    def _help_text(self):
        return """STM Line Profile GUI - Help

PURPOSE
This program loads STMpw-style STM simulation data and creates interactive STM images and line profiles.
It can generate constant-current apparent-height maps, constant-height intensity maps, repeated cells,
non-rectangular unit-cell displays, custom color palettes, and profile data files.

INPUT FILES
Typical input files are:

    TH.dat
    Bardeen.dat
    TH.dat.gz
    Bardeen.dat.gz

The expected input format is:

    N_sampling_z  ngy  ngx
    z             y    x    value
    z             y    x    value
    ...

The first line gives the number of z samples, y grid points, and x grid points. Each following row
contains tip height z, lateral coordinate y, lateral coordinate x, and an STM intensity or
current-like value. The program uses abs(value) internally.

COMPRESSED INPUT FILES
If the file name ends in .gz, the program reads the compressed file directly without creating any
temporary uncompressed file. If you enter TH.dat but only TH.dat.gz exists, the program also detects
and loads the compressed file automatically.

IMAGE MODES
constant current:
    The program searches, at every lateral position, for the height where the requested current setpoint
    is reached. The displayed image is an apparent-height map. Line profiles have apparent height as
    their vertical value.

constant height:
    The program interpolates the STM intensity at the selected height. The displayed image is an intensity
    map. Line profiles have STM intensity as their vertical value.

INTERPOLATION
log:
    Recommended for constant-current STM because tunneling current usually varies approximately
    exponentially with tip height.

linear:
    Direct linear interpolation of the raw current or intensity values.

NUMERIC OPTIONS
Current:
    Current or isosurface setpoint used in constant-current mode.

Height A:
    Tip height in Angstrom used in constant-height mode.

nx, ny:
    Number of periodic repetitions along the first and second lattice-coordinate directions.

UNIT-CELL GEOMETRY
a1_x, a1_y, a2_x, and a2_y define the displayed cell directions. The vectors are normalized internally.
For a rectangular cell use:

    a1_x = 1.0, a1_y = 0.0
    a2_x = 0.0, a2_y = 1.0

For a non-rectangular unit cell, enter the corresponding direction vectors from the original
configuration file. The repeated image is drawn in this transformed geometry, and line-profile
distances are measured in the displayed Cartesian geometry.

PALETTES
The built-in palettes include gray, viridis, plasma, inferno, magma, cividis, hot, afmhot, turbo,
and seismic. Use the Palette dropdown to select a built-in palette or choose load from file... .

load from file...:
    Opens a separate palette browser dialog for gnuplot-style palette files. Clicking a palette file
    previews it and temporarily applies it to the loaded STM image. Press Apply to keep the selected
    custom palette. Press Cancel to restore the previous palette.

LINE PROFILES
Draw a profile:
    Left-click a first point on the STM image, then left-click a second point. The corresponding line
    profile is shown in the profile window. While drawing or dragging a line, a dashed red live-preview
    profile is updated in the profile window. Right-click cancels the current drawing, selection,
    modification, or drag operation. If the profile window is closed accidentally, it is recreated
    automatically the next time profiles are refreshed, for example after adding, selecting, modifying,
    deleting, clearing, undoing, redoing, or restoring lines.

Select line:
    Type 1, 001, or line_001 in the Line box and press Select line. The selected profile statistics
    are shown in the status line and in the profile plot title.

Pick nearest line:
    Press Pick nearest line and click near an existing line on the image. The selected profile
    statistics are shown automatically.

Modify selected:
    Select a line and press Modify selected. The first endpoint is highlighted with a circle and the
    second endpoint with a square. Drag either endpoint to move only that point, drag the line body to
    move the whole line, or click away from the selected line to redefine the line using two clicks.

Delete selected:
    Deletes the selected line.

Clear lines:
    Removes all lines and all profile curves.

SAVING
Save profiles:
    Saves one .dat file per line profile. Constant-current data contain distance and apparent height.
    Constant-height data contain distance and STM intensity.

Save image:
    Saves only the STM image panel as a PNG file, including axes, title, selected/drawn lines, legend,
    and colorbar. The full Matplotlib window can still be saved using the toolbar's built-in save tool.

Save clean:
    Saves the same STM image panel as a PNG file, but without drawn or selected line overlays and
    without the line legend. Axes, title, colorbar, palette, and color range are preserved.

SESSIONS
Session:
    Opens a dialog to save or restore the current working session as a JSON file. A session stores
    the input filename, mode, interpolation, current, height, repetitions, lattice directions, palette,
    selected line, profile counter, and all drawn line coordinates and profile data. Restoring a session
    reloads the image from the original input file and then restores the saved lines.

Ctrl+S:
    Save the current session.

Ctrl+O:
    Restore a saved session.

UNDO AND REDO
Undo / Redo:
    Restore previous or later line-edit states. Undo/redo is available for adding lines, two-click
    modifications, drag modifications, deleting selected lines, clearing lines, and restoring sessions.

Ctrl+Z:
    Undo the last supported edit.

Ctrl+Y or Ctrl+Shift+Z:
    Redo the last undone edit.

PROFILE STATISTICS
When a profile is selected, the program displays the number of valid points, line length, angle,
minimum, maximum, range, mean, and standard deviation. Range is calculated as maximum minus minimum.
The angle is measured in degrees relative to the displayed Cartesian +X direction.

TEXT BOX SHORTCUTS
Ctrl+V or Shift+Insert:
    Paste from clipboard.

Middle mouse button:
    Paste Linux PRIMARY selection, with clipboard fallback.

Tab / Shift+Tab:
    Move to the next or previous text box.

MOUSE CONTROLS AND COORDINATES
Left-click performs the active operation, such as drawing, selecting, modifying, or dragging a line.
Right-click cancels the current operation and returns to draw mode. The status line reports displayed
X and Y coordinates and the image value under the pointer. Long status messages, such as restored
session paths, are shortened in the GUI to avoid covering controls; the full status text is printed to
the terminal.

TROUBLESHOOTING
Blank or missing regions in constant-current mode usually mean the selected current is not reached at
all lateral positions. Try a different current setpoint or inspect the data in constant-height mode.

If constant-height mode fails, check that the selected height lies inside the z range of the input file.

If the displayed cell shape looks wrong, check a1_x, a1_y, a2_x, and a2_y.

If a custom palette does not load, check that the palette file contains a valid set palette defined(...)
block and RGB values between 0 and 1.
"""

    def _show_help(self, event=None):
        help_text = self._help_text()

        if tk is None:
            print(help_text)
            return

        previous_default_root = tk._default_root
        help_root = tk.Tk()
        help_root.title("STM Line Profile GUI Help")
        help_root.geometry("950x720")
        help_root.configure(background="white")

        def close_help():
            try:
                help_root.destroy()
            finally:
                try:
                    if tk._default_root is help_root:
                        tk._default_root = previous_default_root
                except Exception:
                    pass

        help_root.protocol("WM_DELETE_WINDOW", close_help)

        outer = tk.Frame(help_root, background="white")
        outer.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        y_scroll = tk.Scrollbar(outer, orient=tk.VERTICAL)
        y_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        help_textbox = tk.Text(
            outer,
            wrap=tk.WORD,
            yscrollcommand=y_scroll.set,
            background="white",
            foreground="black",
            insertbackground="black",
            selectbackground="#cce8ff",
            selectforeground="black",
            font=("TkFixedFont", 10),
            relief=tk.SOLID,
            borderwidth=1
        )
        help_textbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        y_scroll.config(command=help_textbox.yview)

        help_textbox.insert("1.0", help_text)
        help_textbox.config(state=tk.DISABLED)

        def on_mousewheel(event):
            if getattr(event, "num", None) == 4:
                help_textbox.yview_scroll(-3, "units")
            elif getattr(event, "num", None) == 5:
                help_textbox.yview_scroll(3, "units")
            elif getattr(event, "delta", 0):
                help_textbox.yview_scroll(int(-event.delta / 120), "units")
            return "break"

        help_textbox.bind("<MouseWheel>", on_mousewheel)
        help_textbox.bind("<Button-4>", on_mousewheel)
        help_textbox.bind("<Button-5>", on_mousewheel)
        help_root.bind("<MouseWheel>", on_mousewheel)
        help_root.bind("<Button-4>", on_mousewheel)
        help_root.bind("<Button-5>", on_mousewheel)

        button_frame = tk.Frame(help_root, background="white")
        button_frame.pack(fill=tk.X, padx=8, pady=(0, 8))
        close_button = tk.Button(button_frame, text="Close", command=close_help)
        close_button.pack(side=tk.RIGHT)

        help_root.update_idletasks()
        help_root.lift()
        help_root.focus_force()
        help_root.mainloop()

    def _ensure_colorbar_axes(self):

        if self.ax_cbar is None or self.ax_cbar.figure is None:
            self.ax_cbar = self.fig_img.add_axes(self.cbar_position)
        return self.ax_cbar

    def _clear_colorbar_axes(self):
        if self.cbar is not None:
            try:
                self.cbar.remove()
            except Exception:
                pass
        self.cbar = None
        ax_cbar = self._ensure_colorbar_axes()
        ax_cbar.cla()
        ax_cbar.set_visible(False)

    def _set_mode(self, label):
        self.mode = label

    def _set_interpolation(self, label):
        self.interpolation = label

    def _palette_button_label(self):
        if str(self.palette).startswith("custom:"):
            palette_file = self.custom_palette_file or str(self.palette).split(":", 1)[-1]
            display_name = os.path.basename(palette_file) if palette_file else "custom"
        else:
            display_name = str(self.palette)

        if len(display_name) > 22:
            display_name = display_name[:19] + "..."

        return f"Palette: {display_name}"

    def _update_palette_button_label(self):
        if self.palette_button is None:
            return
        try:
            self.palette_button.label.set_text(self._palette_button_label())
            self.fig_img.canvas.draw_idle()
        except Exception:
            pass

    def _set_palette(self, label):
        if self._ignore_palette_callback:
            return
        if label == "load from file...":
            self._open_custom_palette_dialog()
            return
        self.palette = label
        self.active_cmap = label
        self.custom_palette_file = None
        self.custom_palette_cmap = None
        self._update_palette_button_label()
        self._redraw_palette_now()
        self._set_status(f"palette set to {label}")

    def _redraw_palette_now(self):
        if self.im is None:
            return
        self.im.set_cmap(self.active_cmap)
        if self.cbar is not None:
            try:
                self.cbar.update_normal(self.im)
            except Exception:
                pass
        try:
            self.im.changed()
        except Exception:
            pass
        self.fig_img.canvas.draw()
        self.fig_img.canvas.flush_events()

    def _image_range_text(self):
        if self.image is None:
            return "image range unavailable"

        finite_values = np.asarray(self.image, dtype=float)
        finite_values = finite_values[np.isfinite(finite_values)]

        if finite_values.size == 0:
            return "image range: no finite values"

        image_min = float(np.min(finite_values))
        image_max = float(np.max(finite_values))
        return f"image range: min = {image_min:.6g}, max = {image_max:.6g} {self.image_units}"

    def _open_palette_dropdown_dialog(self, event=None):
        if tk is None or ttk is None:
            self._set_status("tkinter is not available, so the palette dropdown cannot be opened")
            print("tkinter is not available, so the palette dropdown cannot be opened.")
            return

        previous_default_root = tk._default_root
        dialog_root = tk.Tk()
        dialog_root.title("Select palette")
        dialog_root.geometry("390x140")
        dialog_root.configure(background="white")

        selected_palette = tk.StringVar(value=self.palette if self.palette in self.palette_options else "load from file...")

        def close_dialog():
            try:
                dialog_root.destroy()
            finally:
                try:
                    if tk._default_root is dialog_root:
                        tk._default_root = previous_default_root
                except Exception:
                    pass

        selection_applied = {"done": False}

        def apply_selection(event=None):
            if selection_applied["done"]:
                return
            selection_applied["done"] = True
            label = selected_palette.get()
            close_dialog()
            self._set_palette(label)

        dialog_root.protocol("WM_DELETE_WINDOW", close_dialog)

        label_widget = tk.Label(dialog_root, text="Palette:", background="white")
        label_widget.pack(anchor=tk.W, padx=12, pady=(12, 4))

        combo = ttk.Combobox(dialog_root, textvariable=selected_palette, values=self.palette_options, state="readonly", width=34)
        combo.pack(fill=tk.X, padx=12, pady=4)
        combo.bind("<<ComboboxSelected>>", apply_selection)

        button_frame = tk.Frame(dialog_root, background="white")
        button_frame.pack(fill=tk.X, padx=12, pady=10)
        tk.Button(button_frame, text="Apply", command=apply_selection, width=10).pack(side=tk.RIGHT, padx=5)
        tk.Button(button_frame, text="Cancel", command=close_dialog, width=10).pack(side=tk.RIGHT, padx=5)

        dialog_root.update_idletasks()
        dialog_root.lift()
        dialog_root.focus_force()
        combo.focus_set()
        combo.event_generate("<Button-1>")
        dialog_root.mainloop()

    def _get_palette_dropdown_index(self, palette_name):
        for i, name in enumerate(self.palette_options):
            if name == palette_name:
                return i
        return 0

    def _set_palette_dropdown_without_callback(self, palette_name):
        self._ignore_palette_callback = True
        try:
            self.palette = palette_name if palette_name in self.palette_options else "gray"
            self._update_palette_button_label()
        finally:
            self._ignore_palette_callback = False

    def _get_palette_radio_index(self, palette_name):
        return self._get_palette_dropdown_index(palette_name)

    def _set_palette_radio_without_callback(self, palette_name):
        self._set_palette_dropdown_without_callback(palette_name)

    def _open_custom_palette_dialog(self):
        if tk is None or FigureCanvasTkAgg is None:
            self._set_status("tkinter is not available, so palette file dialog cannot be opened")
            return
        previous_palette_name = self.palette
        previous_cmap = self.active_cmap
        initial_directory = os.path.dirname(self.custom_palette_file) if self.custom_palette_file else os.getcwd()
        dialog = CustomPaletteDialog(self, initial_directory, previous_palette_name, previous_cmap)
        dialog.show()

    def _apply_palette_temporarily(self, palette_name, cmap):
        self.palette = palette_name
        self.active_cmap = cmap
        self._update_palette_button_label()
        self._redraw_palette_now()
        self._set_status(f"temporarily previewing palette {palette_name}")

    def _finalize_custom_palette(self, palette_name, palette_file, cmap):
        self.palette = palette_name
        self.active_cmap = cmap
        self.custom_palette_file = palette_file
        self.custom_palette_cmap = cmap
        self._update_palette_button_label()
        self._redraw_palette_now()
        self._set_status(f"custom palette applied: {palette_file}")

    def _restore_palette_after_cancel(self, previous_palette_name, previous_cmap):
        self.palette = previous_palette_name
        self.active_cmap = previous_cmap
        self._update_palette_button_label()
        self._redraw_palette_now()
        if previous_palette_name in self.palette_options:
            self._set_palette_dropdown_without_callback(previous_palette_name)
        else:
            self._set_palette_dropdown_without_callback("gray")
        self._set_status("custom palette selection cancelled. Previous palette restored")

    def _session_key_press(self, event):
        key = str(event.key).lower() if event.key is not None else ""
        if key in ("ctrl+z", "control+z", "cmd+z", "super+z"):
            self._undo()
            return
        if key in ("ctrl+y", "control+y", "cmd+y", "super+y", "ctrl+shift+z", "control+shift+z", "cmd+shift+z", "super+shift+z"):
            self._redo()
            return
        if key in ("ctrl+s", "control+s", "cmd+s", "super+s"):
            self._save_session()
            return
        if key in ("ctrl+o", "control+o", "cmd+o", "super+o"):
            self._restore_session()
            return

    def _open_session_dialog(self, event=None):
        if tk is None:
            print("tkinter is not available. Use Ctrl+S to save to the default session file or Ctrl+O to restore from it.")
            self._set_status("tkinter is not available, session dialog cannot be opened")
            return

        previous_default_root = tk._default_root
        dialog_root = tk.Tk()
        dialog_root.title("Session")
        dialog_root.geometry("360x130")
        dialog_root.configure(background="white")

        def close_dialog():
            try:
                dialog_root.destroy()
            finally:
                try:
                    if tk._default_root is dialog_root:
                        tk._default_root = previous_default_root
                except Exception:
                    pass

        def save_and_close():
            close_dialog()
            self._save_session()

        def restore_and_close():
            close_dialog()
            self._restore_session()

        dialog_root.protocol("WM_DELETE_WINDOW", close_dialog)
        label = tk.Label(dialog_root, text="Save or restore the current STM session.", background="white")
        label.pack(fill=tk.X, padx=10, pady=(12, 8))

        button_frame = tk.Frame(dialog_root, background="white")
        button_frame.pack(fill=tk.X, padx=10, pady=8)
        tk.Button(button_frame, text="Save session", command=save_and_close, width=14).pack(side=tk.LEFT, padx=5)
        tk.Button(button_frame, text="Restore session", command=restore_and_close, width=14).pack(side=tk.LEFT, padx=5)
        tk.Button(button_frame, text="Cancel", command=close_dialog, width=10).pack(side=tk.RIGHT, padx=5)

        dialog_root.update_idletasks()
        dialog_root.lift()
        dialog_root.focus_force()
        dialog_root.mainloop()

    def _default_session_filename(self):
        return f"{self._output_base_name()}_session.json"

    def _choose_session_save_filename(self):
        default_name = self._default_session_filename()
        if tk is None or filedialog is None:
            return default_name

        root = None
        try:
            root = tk.Tk()
            root.withdraw()
            filename = filedialog.asksaveasfilename(
                title="Save STM session",
                initialfile=default_name,
                defaultextension=".json",
                filetypes=(("JSON session files", "*.json"), ("All files", "*.*"))
            )
            return filename if filename else None
        finally:
            if root is not None:
                try:
                    root.destroy()
                except Exception:
                    pass

    def _choose_session_open_filename(self):
        if tk is None or filedialog is None:
            return self._default_session_filename()

        root = None
        try:
            root = tk.Tk()
            root.withdraw()
            filename = filedialog.askopenfilename(
                title="Restore STM session",
                filetypes=(("JSON session files", "*.json"), ("All files", "*.*"))
            )
            return filename if filename else None
        finally:
            if root is not None:
                try:
                    root.destroy()
                except Exception:
                    pass

    def _line_to_session_dict(self, line):
        return {
            "label": line.get("label"),
            "x0": float(line.get("x0", 0.0)),
            "y0": float(line.get("y0", 0.0)),
            "x1": float(line.get("x1", 0.0)),
            "y1": float(line.get("y1", 0.0)),
            "distance": np.asarray(line.get("distance", []), dtype=float).tolist(),
            "profile": np.asarray(line.get("profile", []), dtype=float).tolist(),
            "mode": line.get("mode", self.mode),
            "current": float(line.get("current", self.current_value)),
            "height": float(line.get("height", self.height_value)),
            "units": line.get("units", self.image_units),
            "image_label": line.get("image_label", self.image_label),
            "nx": int(line.get("nx", self.nx)),
            "ny": int(line.get("ny", self.ny)),
            "palette": line.get("palette", self.palette)
        }

    def _line_from_session_dict(self, item):
        return {
            "label": str(item.get("label", "line_000")),
            "x0": float(item.get("x0", 0.0)),
            "y0": float(item.get("y0", 0.0)),
            "x1": float(item.get("x1", 0.0)),
            "y1": float(item.get("y1", 0.0)),
            "distance": np.asarray(item.get("distance", []), dtype=float),
            "profile": np.asarray(item.get("profile", []), dtype=float),
            "mode": str(item.get("mode", self.mode)),
            "current": float(item.get("current", self.current_value)),
            "height": float(item.get("height", self.height_value)),
            "units": str(item.get("units", self.image_units)),
            "image_label": str(item.get("image_label", self.image_label)),
            "nx": int(item.get("nx", self.nx)),
            "ny": int(item.get("ny", self.ny)),
            "palette": str(item.get("palette", self.palette))
        }

    def _make_session_data(self):
        return {
            "session_format": "stm_line_profile_gui_session",
            "session_version": 1,
            "input_file": self.file_box.text.strip(),
            "loaded_source_file": self.data.filename,
            "mode": self.mode,
            "interpolation": self.interpolation,
            "current_value": float(self.current_value),
            "height_value": float(self.height_value),
            "nx": int(self.nx),
            "ny": int(self.ny),
            "a1_x": float(self.a1_x),
            "a1_y": float(self.a1_y),
            "a2_x": float(self.a2_x),
            "a2_y": float(self.a2_y),
            "palette": self.palette,
            "custom_palette_file": self.custom_palette_file,
            "selected_label": self.selected_label,
            "profile_counter": int(self.profile_counter),
            "lines": [self._line_to_session_dict(line) for line in self.lines]
        }

    def _save_session(self, event=None):
        filename = self._choose_session_save_filename()
        if filename is None:
            self._set_status("session save cancelled")
            return

        try:
            session_data = self._make_session_data()
            with open(filename, "w") as f:
                json.dump(session_data, f, indent=2)
            print(f"Saved session {filename}")
            self._set_status(f"session saved: {filename}")
        except Exception as exc:
            print(f"Error while saving session: {exc}")
            self._set_status(f"error saving session: {exc}")

    def _apply_session_controls(self, session_data):
        self.mode = str(session_data.get("mode", self.mode))
        self.interpolation = str(session_data.get("interpolation", self.interpolation))
        self.current_value = float(session_data.get("current_value", self.current_value))
        self.height_value = float(session_data.get("height_value", self.height_value))
        self.nx = int(session_data.get("nx", self.nx))
        self.ny = int(session_data.get("ny", self.ny))
        self.a1_x = float(session_data.get("a1_x", self.a1_x))
        self.a1_y = float(session_data.get("a1_y", self.a1_y))
        self.a2_x = float(session_data.get("a2_x", self.a2_x))
        self.a2_y = float(session_data.get("a2_y", self.a2_y))

        input_file = str(session_data.get("input_file") or session_data.get("loaded_source_file") or self.file_box.text)
        self.file_box.set_val(input_file)
        self.current_box.set_val(f"{self.current_value:.12g}")
        self.height_box.set_val(f"{self.height_value:.12g}")
        self.nx_box.set_val(str(self.nx))
        self.ny_box.set_val(str(self.ny))
        self.a1x_box.set_val(f"{self.a1_x:.12g}")
        self.a1y_box.set_val(f"{self.a1_y:.12g}")
        self.a2x_box.set_val(f"{self.a2_x:.12g}")
        self.a2y_box.set_val(f"{self.a2_y:.12g}")

        if self.mode == "constant current":
            self.mode_radio.set_active(0)
        elif self.mode == "constant height":
            self.mode_radio.set_active(1)

        if self.interpolation == "log":
            self.interp_radio.set_active(0)
        elif self.interpolation == "linear":
            self.interp_radio.set_active(1)

        palette = str(session_data.get("palette", self.palette))
        custom_palette_file = session_data.get("custom_palette_file")
        if palette.startswith("custom:") and custom_palette_file and os.path.exists(custom_palette_file):
            try:
                cmap_name = f"custom_{os.path.splitext(os.path.basename(custom_palette_file))[0]}"
                cmap = load_gnuplot_palette_file(custom_palette_file, cmap_name=cmap_name)
                self._finalize_custom_palette(palette, custom_palette_file, cmap)
            except Exception as exc:
                print(f"Warning: could not restore custom palette {custom_palette_file}: {exc}")
                self._set_palette_radio_without_callback("gray")
                self.palette = "gray"
                self.active_cmap = "gray"
                self.custom_palette_file = None
                self.custom_palette_cmap = None
        elif palette in self.palette_options and palette != "load from file...":
            self._set_palette_radio_without_callback(palette)
            self.palette = palette
            self.active_cmap = palette
            self.custom_palette_file = None
            self.custom_palette_cmap = None
        else:
            self._set_palette_radio_without_callback("gray")
            self.palette = "gray"
            self.active_cmap = "gray"
            self.custom_palette_file = None
            self.custom_palette_cmap = None

        self._update_palette_button_label()

    def _make_edit_state(self, description=""):
        state = self._make_session_data()
        state["history_description"] = description
        state["line_id_text"] = self.line_id_box.text if self.line_id_box is not None else ""
        return state

    def _push_undo_state(self, description):
        if self._restoring_history:
            return

        self.undo_stack.append(self._make_edit_state(description))
        if len(self.undo_stack) > self.undo_limit:
            self.undo_stack = self.undo_stack[-self.undo_limit:]
        self.redo_stack = []

    def _state_matches_current_image(self, state):
        if self.image is None:
            return False

        if state.get("loaded_source_file") != self.data.filename:
            return False

        try:
            checks = [
                self.mode == str(state.get("mode", self.mode)),
                self.interpolation == str(state.get("interpolation", self.interpolation)),
                np.isclose(self.current_value, float(state.get("current_value", self.current_value))),
                np.isclose(self.height_value, float(state.get("height_value", self.height_value))),
                self.nx == int(state.get("nx", self.nx)),
                self.ny == int(state.get("ny", self.ny)),
                np.isclose(self.a1_x, float(state.get("a1_x", self.a1_x))),
                np.isclose(self.a1_y, float(state.get("a1_y", self.a1_y))),
                np.isclose(self.a2_x, float(state.get("a2_x", self.a2_x))),
                np.isclose(self.a2_y, float(state.get("a2_y", self.a2_y)))
            ]
        except Exception:
            return False

        return all(checks)

    def _restore_edit_state(self, state):
        needs_reload = not self._state_matches_current_image(state)
        self._apply_session_controls(state)

        if needs_reload:
            self._load_and_plot()
            if self.image is None:
                return

        self.lines = [self._line_from_session_dict(item) for item in state.get("lines", [])]
        self.profile_counter = int(state.get("profile_counter", len(self.lines)))
        self.selected_label = state.get("selected_label")
        if self.selected_label is not None and self._get_line_by_label(self.selected_label) is None:
            self.selected_label = None

        line_id_text = state.get("line_id_text")
        if line_id_text is None:
            line_id_text = self.selected_label if self.selected_label is not None else ""
        self.line_id_box.set_val(str(line_id_text))

        self.pick_start = None
        self.pending_action = None
        self._clear_live_profile_preview(redraw=False)
        self._clear_drag_state()
        if self.preview_line is not None:
            try:
                self.preview_line.remove()
            except Exception:
                pass
            self.preview_line = None

        self._refresh_all_line_display()

    def _undo(self, event=None):
        if len(self.undo_stack) == 0:
            self._set_status("nothing to undo")
            return

        current_state = self._make_edit_state("redo state")
        state = self.undo_stack.pop()
        self.redo_stack.append(current_state)
        if len(self.redo_stack) > self.undo_limit:
            self.redo_stack = self.redo_stack[-self.undo_limit:]

        self._restoring_history = True
        try:
            self._restore_edit_state(state)
        finally:
            self._restoring_history = False

        description = state.get("history_description", "edit")
        self._set_status(f"undo: {description}")

    def _redo(self, event=None):
        if len(self.redo_stack) == 0:
            self._set_status("nothing to redo")
            return

        current_state = self._make_edit_state("undo state")
        state = self.redo_stack.pop()
        self.undo_stack.append(current_state)
        if len(self.undo_stack) > self.undo_limit:
            self.undo_stack = self.undo_stack[-self.undo_limit:]

        self._restoring_history = True
        try:
            self._restore_edit_state(state)
        finally:
            self._restoring_history = False

        description = state.get("history_description", "edit")
        self._set_status(f"redo: {description}")

    def _restore_session(self, event=None):
        filename = self._choose_session_open_filename()
        if filename is None:
            self._set_status("session restore cancelled")
            return

        try:
            with open(filename, "r") as f:
                session_data = json.load(f)

            if session_data.get("session_format") != "stm_line_profile_gui_session":
                raise ValueError("The selected file is not a valid STM Line Profile GUI session file.")

            if self.image is not None or len(self.lines) > 0:
                self._push_undo_state("restore session")

            self._apply_session_controls(session_data)
            self._load_and_plot()

            restored_lines = [self._line_from_session_dict(item) for item in session_data.get("lines", [])]
            self.lines = restored_lines
            self.profile_counter = int(session_data.get("profile_counter", len(restored_lines)))
            self.selected_label = session_data.get("selected_label")
            if self.selected_label is not None and self._get_line_by_label(self.selected_label) is None:
                self.selected_label = None

            self.line_id_box.set_val(self.selected_label if self.selected_label is not None else "")
            self.pick_start = None
            self.pending_action = None
            if self.preview_line is not None:
                try:
                    self.preview_line.remove()
                except Exception:
                    pass
                self.preview_line = None

            self._refresh_all_line_display()
            image_range_text = self._image_range_text()
            if self.selected_label is not None:
                profile_stats = self._profile_statistics_text(self.selected_label)
                print(profile_stats)
                self._set_status(f"session restored: {filename}. {image_range_text}. {profile_stats}")
            else:
                self._set_status(f"session restored: {filename}. {image_range_text}")
            print(f"Restored session {filename}")
        except Exception as exc:
            print(f"Error while restoring session: {exc}")
            self._set_status(f"error restoring session: {exc}")

    def _prepare_input_file(self, filename):
        filename = filename.strip()

        if filename == "":
            raise FileNotFoundError("No input file was specified.")

        if os.path.exists(filename):
            return filename

        gz_path = filename + ".gz"
        if os.path.exists(gz_path):
            return gz_path

        raise FileNotFoundError(f"File not found: {filename} or {gz_path}")

    def _load_and_plot(self, event=None):
        filename = self.file_box.text.strip()
        try:
            self.current_value = float(self.current_box.text)
            self.height_value = float(self.height_box.text)
            self.nx = int(self.nx_box.text)
            self.ny = int(self.ny_box.text)
            self.a1_x = float(self.a1x_box.text)
            self.a1_y = float(self.a1y_box.text)
            self.a2_x = float(self.a2x_box.text)
            self.a2_y = float(self.a2y_box.text)
        except ValueError:
            print("Could not parse current, height, nx, ny, a1, or a2 values.")
            self._set_status("could not parse current, height, nx, ny, a1, or a2")
            return

        if self.nx < 1 or self.ny < 1:
            print("nx and ny must be positive integers.")
            self._set_status("nx and ny must be positive integers")
            return

        try:
            read_filename = self._prepare_input_file(filename)
            self.data.read(read_filename)
            self.data.filename = read_filename

            if self.mode == "constant current":
                if self.interpolation == "log":
                    self.base_image = constant_current_map_log(self.data.intensity, self.data.z, self.current_value)
                    interp_name = "log interpolation"
                else:
                    self.base_image = constant_current_map_linear(self.data.intensity, self.data.z, self.current_value)
                    interp_name = "linear interpolation"
                self.image_label = f"Constant-current STM, I = {self.current_value:g}, {interp_name}, nx = {self.nx}, ny = {self.ny}"
                self.image_units = "Å"
                cbar_label = "Apparent height (Å)"
            elif self.mode == "constant height":
                self.base_image = constant_height_map(self.data.intensity, self.data.z, self.height_value)
                self.image_label = f"Constant-height STM, z = {self.height_value:g} Å, nx = {self.nx}, ny = {self.ny}"
                self.image_units = "STM intensity"
                cbar_label = "STM intensity"
            else:
                raise ValueError(f"Unknown mode: {self.mode}")

            self.image, self.extended_x, self.extended_y, self.cell_lx, self.cell_ly = extend_periodic_image(
                self.base_image,
                self.data.x,
                self.data.y,
                self.nx,
                self.ny
            )
            self.image_interpolator = RegularGridInterpolator(
                (self.extended_x, self.extended_y),
                self.image,
                bounds_error=False,
                fill_value=np.nan
            )
            self._update_display_geometry()
            self.lines = []
            self.profile_counter = 0
            self.selected_label = None
            self.pick_start = None
            self.pending_action = None
            self.preview_line = None
            self.last_mouse_xy = None
            self._clear_drag_state()
            self.line_id_box.set_val("")
            self._draw_image(cbar_label)
            self._refresh_profile_plot()

            image_range_text = self._image_range_text()
            print(image_range_text)
            self._set_status(f"image loaded. {image_range_text}. Draw, select, modify, or delete lines")

        except Exception as exc:
            print(f"Error while loading or plotting data: {exc}")
            self._set_status(f"error: {exc}")

    def _update_display_geometry(self):
        self.basis, self.inverse_basis = make_lattice_transform(self.a1_x, self.a1_y, self.a2_x, self.a2_y)
        x_edges = axis_edges_from_centers(self.extended_x)
        y_edges = axis_edges_from_centers(self.extended_y)
        u_edges, v_edges = np.meshgrid(x_edges, y_edges, indexing="xy")
        self.cart_x_edges, self.cart_y_edges = lattice_to_cartesian(u_edges, v_edges, self.basis)

    def _cartesian_to_lattice_clamped(self, xcart, ycart):
        u, v = cartesian_to_lattice(xcart, ycart, self.inverse_basis)
        u = np.clip(u, np.min(self.extended_x), np.max(self.extended_x))
        v = np.clip(v, np.min(self.extended_y), np.max(self.extended_y))
        return u, v

    def _lattice_to_cartesian(self, u, v):
        return lattice_to_cartesian(u, v, self.basis)

    def _draw_image(self, cbar_label, update_colorbar=True):
        self.ax_img.clear()
        self.im = self.ax_img.pcolormesh(self.cart_x_edges, self.cart_y_edges, self.image.T, shading="flat", cmap=self.active_cmap)
        self.ax_img.set_xlabel("X (Å)")
        self.ax_img.set_ylabel("Y (Å)")
        self.ax_img.set_aspect("equal")
        self.ax_img.autoscale_view()
        self.ax_img.set_title(self.image_label)

        if update_colorbar or self.cbar is None:
            self._clear_colorbar_axes()
            self.ax_cbar = self._ensure_colorbar_axes()
            self.ax_cbar.set_visible(True)
            self.cbar = self.fig_img.colorbar(self.im, cax=self.ax_cbar)
            self.cbar.set_label(cbar_label)
        else:
            try:
                self.cbar.update_normal(self.im)
                self.cbar.set_label(cbar_label)
            except Exception:
                pass

        self._redraw_image_lines()
        self.fig_img.canvas.draw_idle()

    def _redraw_image_lines(self):
        for line in self.lines:
            label = line["label"]
            is_selected = label == self.selected_label
            color = "red" if is_selected else None
            linewidth = 3.0 if is_selected else 2.0
            self.ax_img.plot([line["x0"], line["x1"]], [line["y0"], line["y1"]], linewidth=linewidth, color=color, label=label)
            if is_selected and self.pending_action == "modify_selected":
                self.ax_img.plot(line["x0"], line["y0"], marker="o", markersize=9, markerfacecolor="white", markeredgecolor="red", markeredgewidth=2.0)
                self.ax_img.plot(line["x1"], line["y1"], marker="s", markersize=9, markerfacecolor="white", markeredgecolor="red", markeredgewidth=2.0)
        if len(self.lines) > 0:
            self.ax_img.legend(loc="upper right", fontsize=8)

    def _profile_statistics(self, label):
        line = self._get_line_by_label(label)
        if line is None:
            return None

        distance, profile = remove_nonfinite_profile_points(line["distance"], line["profile"])
        if len(profile) == 0:
            return {
                "label": label,
                "n": 0,
                "length": 0.0,
                "angle": np.nan,
                "min": np.nan,
                "max": np.nan,
                "range": np.nan,
                "mean": np.nan,
                "std": np.nan,
                "units": line.get("units", self.image_units)
            }

        profile_min = float(np.min(profile))
        profile_max = float(np.max(profile))
        length = float(distance[-1] - distance[0]) if len(distance) > 1 else 0.0
        dx = float(line["x1"] - line["x0"])
        dy = float(line["y1"] - line["y0"])
        angle = float(np.degrees(np.arctan2(dy, dx)))
        return {
            "label": label,
            "n": int(len(profile)),
            "length": length,
            "angle": angle,
            "min": profile_min,
            "max": profile_max,
            "range": profile_max - profile_min,
            "mean": float(np.mean(profile)),
            "std": float(np.std(profile)),
            "units": line.get("units", self.image_units)
        }

    def _profile_statistics_text(self, label):
        stats = self._profile_statistics(label)
        if stats is None:
            return f"{label} stats unavailable"
        if stats["n"] == 0:
            return f"{label} stats: no finite profile values"

        units = stats["units"]
        return (
            f"{label} stats: n = {stats['n']}, "
            f"length = {stats['length']:.6g} Å, "
            f"angle = {stats['angle']:.3f}°, "
            f"min = {stats['min']:.6g}, "
            f"max = {stats['max']:.6g}, "
            f"range = {stats['range']:.6g}, "
            f"mean = {stats['mean']:.6g}, "
            f"std = {stats['std']:.6g} {units}"
        )

    def _show_selected_profile_statistics(self, prefix=None):
        if self.selected_label is None:
            return None

        stats_text = self._profile_statistics_text(self.selected_label)
        print(stats_text)

        if prefix is None or prefix == "":
            self._set_status(stats_text)
        else:
            self._set_status(f"{prefix}. {stats_text}")

        return stats_text

    def _refresh_profile_plot(self):
        self._ensure_profile_window()
        self._clear_live_profile_preview(redraw=False)
        self.ax_prof.clear()
        self.ax_prof.set_xlabel("Distance (Å)")
        self.ax_prof.set_ylabel("Apparent height (Å)" if self.mode == "constant current" else "STM intensity")
        title = "Line profiles"
        if self.selected_label is not None:
            title = f"{title}\n{self._profile_statistics_text(self.selected_label)}"
        self.ax_prof.set_title(title)
        self.ax_prof.grid(True)
        for line in self.lines:
            label = line["label"]
            is_selected = label == self.selected_label
            distance, profile = remove_nonfinite_profile_points(line["distance"], line["profile"])
            if len(distance) < 2:
                continue
            color = "red" if is_selected else None
            linewidth = 2.5 if is_selected else 1.5
            self.ax_prof.plot(distance, profile, marker="o", markersize=3, linewidth=linewidth, color=color, label=label)
        if len(self.lines) > 0:
            self.ax_prof.legend()
        self.fig_prof.canvas.draw_idle()

    def _refresh_image_line_display(self):
        if self.image is not None:
            cbar_label = "Apparent height (Å)" if self.mode == "constant current" else "STM intensity"
            self._draw_image(cbar_label, update_colorbar=False)

    def _refresh_all_line_display(self):
        self._refresh_image_line_display()
        self._refresh_profile_plot()

    def _clamp_to_valid_domain(self, x, y):
        u, v = self._cartesian_to_lattice_clamped(x, y)
        x_clamped, y_clamped = self._lattice_to_cartesian(u, v)
        return float(x_clamped), float(y_clamped)

    def _on_motion(self, event):
        if self.image is None or event.inaxes != self.ax_img or event.xdata is None or event.ydata is None:
            return
        x, y = self._clamp_to_valid_domain(event.xdata, event.ydata)
        self.last_mouse_xy = (x, y)
        value = self._value_at(x, y)
        if np.isfinite(value):
            self.coord_text.set_text(
                f"Mouse: X = {x:.6f} Å, Y = {y:.6f} Å, "
                f"value = {value:.6g} {self.image_units}"
            )
        else:
            self.coord_text.set_text(
                f"Mouse: X = {x:.6f} Å, Y = {y:.6f} Å, value = outside/NaN"
            )
        if self.drag_action is not None:
            self._update_dragged_line_coordinates(x, y)
            line = self._get_line_by_label(self.drag_line_label)
            if line is not None:
                preview_text = self._line_length_angle_text(line["x0"], line["y0"], line["x1"], line["y1"])
                self._update_live_profile_preview(line["x0"], line["y0"], line["x1"], line["y1"])
                if self.drag_action == "start":
                    self._set_status(f"dragging first point of {self.drag_line_label}. {preview_text}")
                elif self.drag_action == "end":
                    self._set_status(f"dragging second point of {self.drag_line_label}. {preview_text}")
                else:
                    self._set_status(f"dragging whole {self.drag_line_label}. {preview_text}")
            self._refresh_image_line_display()
            return

        if self.pick_start is not None:
            x0, y0 = self.pick_start
            if self.preview_line is None:
                line_obj, = self.ax_img.plot([x0, x], [y0, y], linestyle="--", linewidth=1.5, color="red")
                self.preview_line = line_obj
            else:
                self.preview_line.set_data([x0, x], [y0, y])

            preview_text = self._line_length_angle_text(x0, y0, x, y)
            self._update_live_profile_preview(x0, y0, x, y)
            if self.pending_action == "modify_selected":
                self._set_status(f"modifying {self.selected_label}. {preview_text}. Click new second point")
            else:
                self._set_status(f"drawing line. {preview_text}. Click second point")
        self.fig_img.canvas.draw_idle()

    def _get_click_position(self, event):
        if event.inaxes == self.ax_img and event.xdata is not None and event.ydata is not None:
            return self._clamp_to_valid_domain(event.xdata, event.ydata)
        if event.inaxes is not None and event.inaxes != self.ax_img:
            return None
        if self.last_mouse_xy is not None:
            return self.last_mouse_xy
        return None

    def _remove_last_drag_undo_state_if_needed(self):
        if self.drag_line_label is None:
            return

        expected_description = f"drag modify {self.drag_line_label}"
        if len(self.undo_stack) == 0:
            return

        try:
            if self.undo_stack[-1].get("history_description") == expected_description:
                self.undo_stack.pop()
        except Exception:
            pass

    def _cancel_current_operation(self):
        had_operation = False

        if self.drag_action is not None:
            had_operation = True
            self._restore_line_from_drag_copy()
            self._remove_last_drag_undo_state_if_needed()
            self._clear_live_profile_preview(redraw=False)
            self._clear_drag_state()
            self._refresh_all_line_display()
            self._set_status("current operation cancelled. Draw mode")
            return

        if self.pick_start is not None:
            had_operation = True
            self.pick_start = None
            if self.preview_line is not None:
                try:
                    self.preview_line.remove()
                except Exception:
                    pass
                self.preview_line = None

        if self.pending_action is not None:
            had_operation = True
            self.pending_action = None

        if had_operation:
            self._clear_live_profile_preview(redraw=False)
            self._clear_drag_state()
            self._refresh_all_line_display()
            self._set_status("current operation cancelled. Draw mode")
        else:
            self._set_status("no current operation to cancel")

    def _on_click(self, event):
        button = getattr(event, "button", None)

        if button == 3:
            self._cancel_current_operation()
            return

        if button != 1:
            return

        if self.image is None:
            return
        pos = self._get_click_position(event)
        if pos is None:
            return
        x, y = pos
        if self.pending_action == "pick_nearest":
            self._pick_nearest_line_at(x, y)
            self.pending_action = None
            self.pick_start = None
            if self.selected_label is not None:
                self._show_selected_profile_statistics("draw mode")
            else:
                self._set_status("draw mode")
            return
        if self.pending_action == "modify_selected":
            hit_action = self._hit_test_selected_line_for_drag(x, y)
            if hit_action is not None:
                self._start_line_drag(hit_action, x, y)
                return
            self._handle_modify_click(x, y)
            return
        self._handle_draw_click(x, y)

    def _on_release(self, event):
        button = getattr(event, "button", None)

        if button == 3:
            self._cancel_current_operation()
            return

        if button != 1:
            return

        if self.drag_action is None:
            return
        pos = self._get_click_position(event)
        if pos is not None:
            x, y = pos
            self._update_dragged_line_coordinates(x, y)
        self._finalize_dragged_line()

    def _handle_draw_click(self, x, y):
        if self.pick_start is None:
            self.pick_start = (x, y)
            if self.preview_line is not None:
                self.preview_line.remove()
                self.preview_line = None
            self.ax_img.plot(x, y, marker="o", markersize=5, color="red")
            self._clear_live_profile_preview()
            self._ensure_profile_window()
            self.fig_img.canvas.draw_idle()
            preview_text = self._line_length_angle_text(x, y, x, y)
            self._set_status(f"first point selected. {preview_text}. Click second point")
            return
        x0, y0 = self.pick_start
        x1, y1 = x, y
        if self.preview_line is not None:
            self.preview_line.remove()
            self.preview_line = None
        distance, profile = self._extract_profile(x0, y0, x1, y1)
        if len(distance) < 2:
            print("Selected line has fewer than two valid finite profile points. Line was not added.")
            self.pick_start = None
            self._set_status("line not added because it has too few valid points")
            self.fig_img.canvas.draw_idle()
            return
        self._push_undo_state("add line")
        self.profile_counter += 1
        label = f"line_{self.profile_counter:03d}"
        self.lines.append({"label": label, "x0": x0, "y0": y0, "x1": x1, "y1": y1, "distance": distance, "profile": profile, "mode": self.mode, "current": self.current_value, "height": self.height_value, "units": self.image_units, "image_label": self.image_label, "nx": self.nx, "ny": self.ny, "palette": self.palette})
        self.selected_label = label
        self.line_id_box.set_val(label)
        self.pick_start = None
        self._clear_live_profile_preview(redraw=False)
        self._refresh_all_line_display()
        self._show_selected_profile_statistics(f"added and selected {label}")

    def _handle_modify_click(self, x, y):
        if self.selected_label is None:
            self.pending_action = None
            self.pick_start = None
            self._set_status("no selected line to modify")
            return
        if self.pick_start is None:
            self.pick_start = (x, y)
            if self.preview_line is not None:
                self.preview_line.remove()
                self.preview_line = None
            self.ax_img.plot(x, y, marker="o", markersize=5, color="red")
            self._clear_live_profile_preview()
            self._ensure_profile_window()
            self.fig_img.canvas.draw_idle()
            preview_text = self._line_length_angle_text(x, y, x, y)
            self._set_status(f"modifying {self.selected_label}. {preview_text}. Click new second point")
            return
        x0, y0 = self.pick_start
        x1, y1 = x, y
        if self.preview_line is not None:
            self.preview_line.remove()
            self.preview_line = None
        distance, profile = self._extract_profile(x0, y0, x1, y1)
        if len(distance) < 2:
            print("Modified line has fewer than two valid finite profile points. Modification was cancelled.")
            self.pick_start = None
            self.pending_action = None
            self._set_status("modification cancelled because line has too few valid points")
            self._refresh_all_line_display()
            return
        line = self._get_line_by_label(self.selected_label)
        if line is None:
            self.pick_start = None
            self.pending_action = None
            self._set_status("selected line was not found")
            self._refresh_all_line_display()
            return
        self._push_undo_state(f"modify {self.selected_label}")
        line.update({"x0": x0, "y0": y0, "x1": x1, "y1": y1, "distance": distance, "profile": profile, "mode": self.mode, "current": self.current_value, "height": self.height_value, "units": self.image_units, "image_label": self.image_label, "nx": self.nx, "ny": self.ny, "palette": self.palette})
        self.pick_start = None
        self.pending_action = None
        self._clear_live_profile_preview(redraw=False)
        self._refresh_all_line_display()
        self._show_selected_profile_statistics(f"modified {self.selected_label}")

    def _extract_profile(self, x0, y0, x1, y1):
        x0, y0 = self._clamp_to_valid_domain(x0, y0)
        x1, y1 = self._clamp_to_valid_domain(x1, y1)
        dx = x1 - x0
        dy = y1 - y0
        length = np.sqrt(dx * dx + dy * dy)
        if length == 0.0:
            return np.array([]), np.array([])
        step_u = abs(np.mean(np.diff(self.extended_x)))
        step_v = abs(np.mean(np.diff(self.extended_y)))
        nominal_step = min(step_u, step_v)
        npts = max(2, int(np.ceil(length / nominal_step)) + 1)
        xcart = np.linspace(x0, x1, npts)
        ycart = np.linspace(y0, y1, npts)
        u, v = cartesian_to_lattice(xcart, ycart, self.inverse_basis)
        u = np.clip(u, np.min(self.extended_x), np.max(self.extended_x))
        v = np.clip(v, np.min(self.extended_y), np.max(self.extended_y))
        xcart, ycart = self._lattice_to_cartesian(u, v)
        profile = self.image_interpolator(np.column_stack((u, v)))
        distance = np.sqrt((xcart - xcart[0]) ** 2 + (ycart - ycart[0]) ** 2)
        return remove_nonfinite_profile_points(distance, profile)

    def _value_at(self, x, y):
        if self.image_interpolator is None:
            return np.nan
        u, v = self._cartesian_to_lattice_clamped(x, y)
        return self.image_interpolator([[u, v]])[0]

    def _normalize_line_label(self, text):
        text = text.strip()
        if text == "":
            return None
        if text.startswith("line_"):
            return text
        try:
            number = int(text)
            return f"line_{number:03d}"
        except ValueError:
            return text

    def _get_line_by_label(self, label):
        for line in self.lines:
            if line["label"] == label:
                return line
        return None

    def _select_line_from_textbox(self, event=None):
        label = self._normalize_line_label(self.line_id_box.text)
        if label is None:
            self._set_status("no line id entered")
            return
        line = self._get_line_by_label(label)
        if line is None:
            self._set_status(f"{label} not found")
            print(f"{label} not found.")
            return
        self.selected_label = label
        self.line_id_box.set_val(label)
        self._refresh_all_line_display()
        self._show_selected_profile_statistics(f"selected {label}")

    def _start_pick_nearest_line(self, event=None):
        if len(self.lines) == 0:
            self._set_status("no lines available to select")
            print("No lines available to select.")
            return
        self.pending_action = "pick_nearest"
        self.pick_start = None
        self._set_status("click near a line in the STM image")

    def _pick_nearest_line_at(self, x, y):
        if len(self.lines) == 0:
            self._set_status("no lines available")
            return
        xmin = float(np.nanmin(self.cart_x_edges))
        xmax = float(np.nanmax(self.cart_x_edges))
        ymin = float(np.nanmin(self.cart_y_edges))
        ymax = float(np.nanmax(self.cart_y_edges))
        diag = np.sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2)
        step_x = abs(np.mean(np.diff(self.extended_x)))
        step_y = abs(np.mean(np.diff(self.extended_y)))
        tolerance = max(3.0 * min(step_x, step_y), 0.01 * diag)
        best_label = None
        best_distance = None
        for line in self.lines:
            d = point_to_segment_distance(x, y, line["x0"], line["y0"], line["x1"], line["y1"])
            if best_distance is None or d < best_distance:
                best_distance = d
                best_label = line["label"]
        if best_label is None:
            self._set_status("no line found")
            return
        if best_distance is not None and best_distance > tolerance:
            print(f"Nearest line is {best_label}, but distance {best_distance:.4f} Å is larger than selection tolerance {tolerance:.4f} Å.")
            self._set_status("click closer to a line")
            return
        self.selected_label = best_label
        self.line_id_box.set_val(best_label)
        self._refresh_all_line_display()
        self._show_selected_profile_statistics(f"selected {best_label}")

    def _start_modify_selected_line(self, event=None):
        if self.selected_label is None:
            self._set_status("select a line before modifying")
            print("Select a line before modifying.")
            return
        if self._get_line_by_label(self.selected_label) is None:
            self._set_status("selected line not found")
            print("Selected line not found.")
            return
        self.pending_action = "modify_selected"
        self.pick_start = None
        self._clear_drag_state()
        self._refresh_all_line_display()
        self._set_status(f"modifying {self.selected_label}. Drag an endpoint, drag the line body, or click away to redefine the line")

    def _delete_selected_line(self, event=None):
        label = self.selected_label if self.selected_label is not None else self._normalize_line_label(self.line_id_box.text)
        if label is None:
            self._set_status("no selected line to delete")
            print("No selected line to delete.")
            return
        old_count = len(self.lines)
        self._push_undo_state(f"delete {label}")
        self.lines = [line for line in self.lines if line["label"] != label]
        if len(self.lines) == old_count:
            self._set_status(f"{label} not found")
            print(f"{label} not found.")
            return
        self.selected_label = None
        self.line_id_box.set_val("")
        self.pick_start = None
        self.pending_action = None
        self._clear_live_profile_preview(redraw=False)
        self._clear_drag_state()
        self._refresh_all_line_display()
        self._set_status(f"deleted {label}")

    def _clear_lines(self, event=None):
        if len(self.lines) == 0:
            self._set_status("no lines to clear")
            return
        self._push_undo_state("clear lines")
        self.lines = []
        self.profile_counter = 0
        self.selected_label = None
        self.pick_start = None
        self.pending_action = None
        self._clear_drag_state()
        if self.preview_line is not None:
            self.preview_line.remove()
            self.preview_line = None
        self.line_id_box.set_val("")
        self._refresh_all_line_display()
        self._set_status("all lines cleared. Draw mode")

    def _save_profiles(self, event=None):
        if len(self.lines) == 0:
            print("No line profiles to save.")
            self._set_status("no line profiles to save")
            return
        base = self._output_base_name()
        for line in self.lines:
            distance, profile = remove_nonfinite_profile_points(line["distance"], line["profile"])
            if len(distance) == 0:
                print(f"Profile {line['label']} has no valid finite points. File was not saved.")
                continue
            filename = f"{base}_{line['label']}.dat"
            with open(filename, "w") as f:
                f.write("# STM line profile generated from STMpw data\n")
                f.write(f"# Source file: {self.data.filename}\n")
                f.write(f"# Image mode: {line['mode']}\n")
                f.write(f"# Lateral repetition: nx = {line['nx']}, ny = {line['ny']}\n")
                f.write(f"# Lattice directions: a1 = ({self.a1_x:.12g}, {self.a1_y:.12g}), a2 = ({self.a2_x:.12g}, {self.a2_y:.12g})\n")
                f.write("# Line coordinates are Cartesian display coordinates after applying the lattice directions.\n")
                f.write(f"# Palette used for image display: {line['palette']}\n")
                if line["mode"] == "constant current":
                    f.write(f"# Current setpoint: {line['current']:.12g}\n")
                    f.write("# Profile value: apparent height in Angstrom\n")
                    f.write("# Columns: distance_Angstrom apparent_height_Angstrom\n")
                else:
                    f.write(f"# Constant height: {line['height']:.12g} Angstrom\n")
                    f.write("# Profile value: STM intensity in native STMpw units\n")
                    f.write("# Columns: distance_Angstrom STM_intensity\n")
                f.write(f"# Line start: x0 = {line['x0']:.12g} Angstrom, y0 = {line['y0']:.12g} Angstrom\n")
                f.write(f"# Line end:   x1 = {line['x1']:.12g} Angstrom, y1 = {line['y1']:.12g} Angstrom\n")
                f.write("# Non-finite values were removed before saving.\n")
                for d, val in zip(distance, profile):
                    f.write(f"{d:20.10f} {val:20.10e}\n")
            print(f"Saved {filename}")
        self._set_status("profiles saved")

    def _save_figures(self, event=None):
        if self.image is None or self.im is None:
            print("No STM image to save.")
            self._set_status("no STM image to save")
            return

        base = self._output_base_name()
        image_png = f"{base}_stm_image.png"

        try:
            self.fig_img.canvas.draw()
            renderer = self.fig_img.canvas.get_renderer()
            bboxes = [self.ax_img.get_tightbbox(renderer)]

            if self.ax_cbar is not None and self.ax_cbar.get_visible():
                bboxes.append(self.ax_cbar.get_tightbbox(renderer))

            bbox = Bbox.union(bboxes)
            bbox = bbox.transformed(self.fig_img.dpi_scale_trans.inverted())
            bbox = bbox.expanded(1.04, 1.08)
            self.fig_img.savefig(image_png, dpi=300, bbox_inches=bbox)
        except Exception as exc:
            print(f"Warning: precise STM image crop failed ({exc}). Saving image axes fallback.")
            extent = self.ax_img.get_window_extent().transformed(self.fig_img.dpi_scale_trans.inverted())
            extent = extent.expanded(1.08, 1.12)
            self.fig_img.savefig(image_png, dpi=300, bbox_inches=extent)

        print(f"Saved {image_png}")
        self._set_status(f"STM image saved: {image_png}")

    def _save_clean_image(self, event=None):
        if self.image is None:
            print("No STM image to save.")
            self._set_status("no STM image to save")
            return

        base = self._output_base_name()
        image_png = f"{base}_stm_image_clean.png"
        cbar_label = "Apparent height (Å)" if self.mode == "constant current" else "STM intensity"

        fig_clean, ax_clean = plt.subplots(figsize=(8.5, 6.0))
        try:
            clean_im = ax_clean.pcolormesh(
                self.cart_x_edges,
                self.cart_y_edges,
                self.image.T,
                shading="flat",
                cmap=self.active_cmap
            )

            if self.im is not None:
                try:
                    clean_im.set_clim(*self.im.get_clim())
                except Exception:
                    pass

            ax_clean.set_xlabel("X (Å)")
            ax_clean.set_ylabel("Y (Å)")
            ax_clean.set_aspect("equal")
            ax_clean.autoscale_view()
            ax_clean.set_title(self.image_label)
            cbar = fig_clean.colorbar(clean_im, ax=ax_clean)
            cbar.set_label(cbar_label)
            fig_clean.savefig(image_png, dpi=300, bbox_inches="tight")
        finally:
            plt.close(fig_clean)

        print(f"Saved {image_png}")
        self._set_status(f"clean STM image saved: {image_png}")

    def _output_base_name(self):
        source = self.data.filename if self.data.filename is not None else "stm"
        basename = os.path.basename(source)
        if basename.endswith(".gz"):
            basename = basename[:-3]
        root = os.path.splitext(basename)[0]
        if self.mode == "constant current":
            tag = f"cc_I_{self.current_value:.6g}".replace("-", "m").replace("+", "")
        else:
            tag = f"ch_z_{self.height_value:.6g}".replace("-", "m").replace("+", "")
        tag = tag.replace(".", "p").replace("e", "E")
        return f"{root}_{tag}_nx{self.nx}_ny{self.ny}"


def main():
    gui = STMLineProfileGUI()
    plt.show()


if __name__ == "__main__":
    main()
