/****************************************************************************
** Meta object code from reading C++ file 'UfabMainWindow.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.7.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "UfabMainWindow.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'UfabMainWindow.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.7.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_UfabMainWindow_t {
    QByteArrayData data[7];
    char stringdata0[64];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_UfabMainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_UfabMainWindow_t qt_meta_stringdata_UfabMainWindow = {
    {
        QT_MOC_LITERAL(0, 0, 14), // "UfabMainWindow"
        QT_MOC_LITERAL(1, 15, 14), // "openFileDialog"
        QT_MOC_LITERAL(2, 30, 0), // ""
        QT_MOC_LITERAL(3, 31, 11), // "resizeEvent"
        QT_MOC_LITERAL(4, 43, 13), // "QResizeEvent*"
        QT_MOC_LITERAL(5, 57, 1), // "e"
        QT_MOC_LITERAL(6, 59, 4) // "Quit"

    },
    "UfabMainWindow\0openFileDialog\0\0"
    "resizeEvent\0QResizeEvent*\0e\0Quit"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_UfabMainWindow[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    3,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   29,    2, 0x08 /* Private */,
    3,    1,   30,    2, 0x08 /* Private */,
    6,    0,   33,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 4,    5,
    QMetaType::Void,

    0        // eod
};

void UfabMainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        UfabMainWindow *_t = static_cast<UfabMainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->openFileDialog();
            break;
        case 1:
            _t->resizeEvent((*reinterpret_cast< QResizeEvent*(*)>(_a[1])));
            break;
        case 2:
            _t->Quit();
            break;
        default:
            ;
        }
    }
}

const QMetaObject UfabMainWindow::staticMetaObject = {
    {   &QMainWindow::staticMetaObject, qt_meta_stringdata_UfabMainWindow.data,
        qt_meta_data_UfabMainWindow,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *UfabMainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *UfabMainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_UfabMainWindow.stringdata0))
        return static_cast<void*>(const_cast< UfabMainWindow*>(this));
    if (!strcmp(_clname, "Ui::UfabMainWindow"))
        return static_cast< Ui::UfabMainWindow*>(const_cast< UfabMainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int UfabMainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
