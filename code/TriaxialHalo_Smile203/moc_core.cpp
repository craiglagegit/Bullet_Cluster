/****************************************************************************
** Meta object code from reading C++ file 'core.h'
**
** Created: Sat Aug 25 08:31:27 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "core.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'core.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_smile__CCalcThread[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      20,   19,   19,   19, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_smile__CCalcThread[] = {
    "smile::CCalcThread\0\0calcFinished()\0"
};

const QMetaObject smile::CCalcThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_smile__CCalcThread,
      qt_meta_data_smile__CCalcThread, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &smile::CCalcThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *smile::CCalcThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *smile::CCalcThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_smile__CCalcThread))
        return static_cast<void*>(const_cast< CCalcThread*>(this));
    return QThread::qt_metacast(_clname);
}

int smile::CCalcThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: calcFinished(); break;
        default: ;
        }
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void smile::CCalcThread::calcFinished()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
static const uint qt_meta_data_smile__CCalcManyThread[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      24,   23,   23,   23, 0x05,

 // slots: signature, parameters, type, tag, flags
      39,   23,   23,   23, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_smile__CCalcManyThread[] = {
    "smile::CCalcManyThread\0\0calcFinished()\0"
    "stopThread()\0"
};

const QMetaObject smile::CCalcManyThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_smile__CCalcManyThread,
      qt_meta_data_smile__CCalcManyThread, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &smile::CCalcManyThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *smile::CCalcManyThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *smile::CCalcManyThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_smile__CCalcManyThread))
        return static_cast<void*>(const_cast< CCalcManyThread*>(this));
    return QThread::qt_metacast(_clname);
}

int smile::CCalcManyThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: calcFinished(); break;
        case 1: stopThread(); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void smile::CCalcManyThread::calcFinished()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
static const uint qt_meta_data_smile__CSchwarzschildThread[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      36,   29,   28,   28, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_smile__CSchwarzschildThread[] = {
    "smile::CSchwarzschildThread\0\0result\0"
    "calcFinished(QString)\0"
};

const QMetaObject smile::CSchwarzschildThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_smile__CSchwarzschildThread,
      qt_meta_data_smile__CSchwarzschildThread, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &smile::CSchwarzschildThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *smile::CSchwarzschildThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *smile::CSchwarzschildThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_smile__CSchwarzschildThread))
        return static_cast<void*>(const_cast< CSchwarzschildThread*>(this));
    return QThread::qt_metacast(_clname);
}

int smile::CSchwarzschildThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: calcFinished((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void smile::CSchwarzschildThread::calcFinished(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_smile__CNbodyExportStatThread[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      38,   31,   30,   30, 0x05,

 // slots: signature, parameters, type, tag, flags
      62,   30,   30,   30, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_smile__CNbodyExportStatThread[] = {
    "smile::CNbodyExportStatThread\0\0result\0"
    "exportFinished(QString)\0stopThread()\0"
};

const QMetaObject smile::CNbodyExportStatThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_smile__CNbodyExportStatThread,
      qt_meta_data_smile__CNbodyExportStatThread, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &smile::CNbodyExportStatThread::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *smile::CNbodyExportStatThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *smile::CNbodyExportStatThread::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_smile__CNbodyExportStatThread))
        return static_cast<void*>(const_cast< CNbodyExportStatThread*>(this));
    return QThread::qt_metacast(_clname);
}

int smile::CNbodyExportStatThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: exportFinished((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: stopThread(); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void smile::CNbodyExportStatThread::exportFinished(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_smile__CSmileCore[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
      20,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       9,       // signalCount

 // signals: signature, parameters, type, tag, flags
      24,   19,   18,   18, 0x05,
      44,   18,   18,   18, 0x05,
      58,   18,   18,   18, 0x05,
      80,   18,   18,   18, 0x05,
     104,   18,   18,   18, 0x05,
     144,  137,   18,   18, 0x05,
     184,  137,   18,   18, 0x05,
     219,   18,   18,   18, 0x05,
     236,   18,   18,   18, 0x05,

 // slots: signature, parameters, type, tag, flags
     249,   18,   18,   18, 0x08,
     269,   18,   18,   18, 0x08,
     291,   18,   18,   18, 0x08,
     322,  137,   18,   18, 0x08,
     360,   18,   18,   18, 0x08,
     378,   19,   18,   18, 0x08,
     411,   18,   18,   18, 0x08,
     431,   18,   18,   18, 0x08,
     468,  460,   18,   18, 0x08,
     504,   18,   18,   18, 0x08,
     523,  460,   18,   18, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_smile__CSmileCore[] = {
    "smile::CSmileCore\0\0info\0signalInfo(QString)\0"
    "signalTimer()\0signalOrbitFinished()\0"
    "signalFreqMapFinished()\0"
    "signalSchwOrbitLibraryFinished()\0"
    "result\0signalSchwOptimizationFinished(QString)\0"
    "signalNbodyExportFinished(QString)\0"
    "scriptNextLine()\0scriptDone()\0"
    "coreOrbitFinished()\0coreFreqMapFinished()\0"
    "coreSchwOrbitLibraryFinished()\0"
    "coreSchwOptimizationFinished(QString)\0"
    "coreNbodyExport()\0coreNbodyExportFinished(QString)\0"
    "scriptProcessLine()\0scriptOrbitLibraryFinished()\0"
    "message\0scriptOptimizationFinished(QString)\0"
    "scriptTimerEvent()\0scriptInfo(QString)\0"
};

const QMetaObject smile::CSmileCore::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_smile__CSmileCore,
      qt_meta_data_smile__CSmileCore, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &smile::CSmileCore::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *smile::CSmileCore::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *smile::CSmileCore::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_smile__CSmileCore))
        return static_cast<void*>(const_cast< CSmileCore*>(this));
    return QObject::qt_metacast(_clname);
}

int smile::CSmileCore::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: signalInfo((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: signalTimer(); break;
        case 2: signalOrbitFinished(); break;
        case 3: signalFreqMapFinished(); break;
        case 4: signalSchwOrbitLibraryFinished(); break;
        case 5: signalSchwOptimizationFinished((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: signalNbodyExportFinished((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 7: scriptNextLine(); break;
        case 8: scriptDone(); break;
        case 9: coreOrbitFinished(); break;
        case 10: coreFreqMapFinished(); break;
        case 11: coreSchwOrbitLibraryFinished(); break;
        case 12: coreSchwOptimizationFinished((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 13: coreNbodyExport(); break;
        case 14: coreNbodyExportFinished((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 15: scriptProcessLine(); break;
        case 16: scriptOrbitLibraryFinished(); break;
        case 17: scriptOptimizationFinished((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 18: scriptTimerEvent(); break;
        case 19: scriptInfo((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 20;
    }
    return _id;
}

// SIGNAL 0
void smile::CSmileCore::signalInfo(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void smile::CSmileCore::signalTimer()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void smile::CSmileCore::signalOrbitFinished()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}

// SIGNAL 3
void smile::CSmileCore::signalFreqMapFinished()
{
    QMetaObject::activate(this, &staticMetaObject, 3, 0);
}

// SIGNAL 4
void smile::CSmileCore::signalSchwOrbitLibraryFinished()
{
    QMetaObject::activate(this, &staticMetaObject, 4, 0);
}

// SIGNAL 5
void smile::CSmileCore::signalSchwOptimizationFinished(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void smile::CSmileCore::signalNbodyExportFinished(const QString & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void smile::CSmileCore::scriptNextLine()
{
    QMetaObject::activate(this, &staticMetaObject, 7, 0);
}

// SIGNAL 8
void smile::CSmileCore::scriptDone()
{
    QMetaObject::activate(this, &staticMetaObject, 8, 0);
}
QT_END_MOC_NAMESPACE
