# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Geo2Local
                                 A QGIS plugin
 This plugin converts coordinates
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2018-04-19
        git sha              : $Format:%H$
        copyright            : (C) 2018 by Diego Moreira Carvalho/Curupira Tecnologia
        email                : diego@curupiratecnologia.com.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt5.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QAction, QFileDialog, QMessageBox

# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .geo2local_dialog import Geo2LocalDialog
import os.path
from .converter.ellipsoid import ellipsoids, elipsoidsIndex, Ellipsoid
from .converter.transformation import TransformationRT, TransformationNBR14166

class Geo2Local:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'Geo2Local_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = Geo2LocalDialog()
        for k,v in elipsoidsIndex.items():
            self.dlg.elipsoideComboBox.addItem(ellipsoids[v]['description'])

        
        self.dlg.entradaPushButton.clicked.connect(self.selecionarEntrada)

        self.dlg.saidaPushButton.clicked.connect(self.selecionarSaida)

        self.dlg.button_box.accepted.connect(self.processar)


        self.dlg.longLineEdit.setText('-53.8033814')
        self.dlg.latLineEdit.setText('-29.6851191')
        self.dlg.altLineEdit.setText('135.788')

        self.dlg.xLineEdit.setText('150000')
        self.dlg.yLineEdit.setText('250000')
        
        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Geo2Local')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'Geo2Local')
        self.toolbar.setObjectName(u'Geo2Local')


        
    def processar(self):
        fileIn = self.dlg.entradaLineEdit.text()
        fileOut = self.dlg.saidaLineEdit.text()

        long = float(self.dlg.longLineEdit.text())
        lat = float(self.dlg.latLineEdit.text())
        alt = float(self.dlg.altLineEdit.text())

        dX = float(self.dlg.xLineEdit.text())
        dY = float(self.dlg.yLineEdit.text())      
        
        

        elipsoide = self.dlg.elipsoideComboBox.currentIndex() + 1
        metodo = self.dlg.metodoComboBox.currentIndex()
        direcao = self.dlg.direcaoComboBox.currentIndex()

        origin = [lat, long, alt]

        if metodo == 0 and direcao == 0:
            QMessageBox.warning(self.dlg, "Aviso","Ainda não Implementado!!!")
            
        else:
            el = Ellipsoid(name=elipsoidsIndex[elipsoide])
            if metodo == 0:
                t = TransformationNBR14166(el, origin)
            else:
                t = TransformationRT(el, origin)


            entrada = open(fileIn)
            lines = [i for i in entrada.readlines() if 'nome' not in i]
            if direcao == 0:
                coords = [t.local2geo([float(j) for j in i.replace('\n','').split(',')[1:]]) for i in lines]
            else:
                coords = [t.geo2local([float(j) for j in i.replace('\n','').split(',')[1:]]) for i in lines]
                coords = [[i[0]+dY, i[1]+dX, i[2]] for i in coords]

            indice = [i.replace('\n','').split(',')[0] for i in lines]

            saida = open(fileOut, 'w')
            saida.writelines('nome,lat,long,alt\n')
            for i in range(len(indice)):
                line = '%s,%s\n' % (indice[i], ','.join([str(i) for i in coords[i]]))
                saida.writelines(line)
            entrada.close()
            saida.close()


    def selecionarEntrada(self):
        file, _ = QFileDialog.getOpenFileName(self.dlg, "Selecione um arquivo", filter='*.csv')
        self.dlg.entradaLineEdit.setText(file)


    def selecionarSaida(self):
        file, _ = QFileDialog.getSaveFileName(self.dlg, "Selecione um arquivo", filter='*.csv')
        self.dlg.saidaLineEdit.setText(file)


    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('Geo2Local', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToVectorMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/geo2local/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Geo2Local'),
            callback=self.run,
            parent=self.iface.mainWindow())


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginVectorMenu(
                self.tr(u'&Geo2Local'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar


    def run(self):
        """Run method that performs all the real work"""
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            # Do something useful here - delete the line containing pass and
            # substitute with your code.
            pass
